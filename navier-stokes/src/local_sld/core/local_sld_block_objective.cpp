#include "local_sld_block_objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

TriadSelection complementary_selection(TriadSelection selection) {
    using Mode = TriadSelection::SignatureMode;
    switch (selection.signature_mode()) {
        case Mode::include: {
            const auto signature = selection.squared_length_signature();
            return TriadSelection::local_without_signature(
                signature[0], signature[1], signature[2]);
        }
        case Mode::exclude: {
            const auto signature = selection.squared_length_signature();
            return TriadSelection::local_signature(
                signature[0], signature[1], signature[2]);
        }
        case Mode::include_equal_low_doubling:
            return TriadSelection::local_without_equal_low_doubling();
        case Mode::exclude_equal_low_doubling:
            return TriadSelection::local_equal_low_doubling();
        case Mode::exclude_equal_low_doubling_and_signature:
            break;
        case Mode::include_equal_low_double_triple:
            return TriadSelection::local_without_equal_low_double_triple();
        case Mode::exclude_equal_low_double_triple:
            return TriadSelection::local_equal_low_double_triple();
        case Mode::exclude_equal_low_double_triple_and_signature:
            break;
        case Mode::none:
            break;
    }
    throw std::invalid_argument(
        "SLD block objective requires a signature selection or its complement");
}

void scale(SpectralIncrement& value, SpectralReal factor) {
    for (ComplexVector& mode : value) {
        for (SpectralComplex& component : mode) {
            component *= factor;
        }
    }
}

void add_scaled(SpectralIncrement& target,
                const SpectralIncrement& source,
                SpectralReal factor) {
    if (target.size() != source.size()) {
        throw std::invalid_argument("SLD block gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += factor * source[mode][component];
        }
    }
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

SpectralReal relative_error(SpectralReal left, SpectralReal right) {
    return std::abs(left - right) /
        std::max({std::abs(left), std::abs(right), 1e-30L});
}

SpectralReal select_constant(
    const LocalSldBlockObjectiveValue& value,
    LocalSldBlock block) {
    switch (block) {
        case LocalSldBlock::full:
            return value.full_constant_ratio;
        case LocalSldBlock::selected_closed:
            return value.selected_constant_ratio;
        case LocalSldBlock::complement_closed:
            return value.complement_constant_ratio;
        case LocalSldBlock::mixed:
            return value.mixed_constant_ratio;
    }
    return 0.0L;
}

}  // namespace

LocalSldBlockObjective::LocalSldBlockObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selected,
    LocalSldBlock block)
    : dynamics_(dynamics), selected_(selected),
      complement_(complementary_selection(selected)), block_(block) {}

LocalSldBlockObjectiveValue LocalSldBlockObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjective full_objective(
        dynamics_, TriadPartition::local);
    const LocalQuarticClosureObjective selected_objective(
        dynamics_, selected_);
    const LocalQuarticClosureObjective complement_objective(
        dynamics_, complement_);
    const LocalQuarticClosureObjectiveValue full =
        full_objective.evaluate(state);
    const LocalQuarticClosureObjectiveValue selected =
        selected_objective.evaluate(state);
    const LocalQuarticClosureObjectiveValue complement =
        complement_objective.evaluate(state);

    LocalSldBlockObjectiveValue result;
    result.full_constant_ratio = full.signed_constant_ratio;
    result.selected_constant_ratio = selected.signed_constant_ratio;
    result.complement_constant_ratio = complement.signed_constant_ratio;
    result.mixed_constant_ratio = result.full_constant_ratio -
        result.selected_constant_ratio -
        result.complement_constant_ratio;
    result.block_constant_ratio = select_constant(result, block_);
    result.normalized_stretching =
        full.normalized_stretching_ratio;
    result.common_shape_factor = full.signed_shape_factor;
    result.block_sld_ratio = result.block_constant_ratio *
        result.common_shape_factor;
    result.full_sld_ratio = full.signed_local_sld_ratio;
    const SpectralReal reconstructed =
        (result.selected_constant_ratio +
         result.complement_constant_ratio +
         result.mixed_constant_ratio) *
        result.common_shape_factor;
    result.ratio_reconstruction_error = relative_error(
        result.full_sld_ratio, reconstructed);
    result.finite = full.finite && selected.finite && complement.finite &&
        std::isfinite(result.block_sld_ratio) &&
        result.ratio_reconstruction_error < 1e-12L;
    return result;
}

SpectralIncrement LocalSldBlockObjective::gradient(
    const SpectralState& state) const {
    const LocalQuarticClosureObjective full_objective(
        dynamics_, TriadPartition::local);
    const LocalQuarticClosureObjective selected_objective(
        dynamics_, selected_);
    const LocalQuarticClosureObjective complement_objective(
        dynamics_, complement_);
    const LocalQuarticClosureObjectiveValue full =
        full_objective.evaluate(state);
    const LocalSldBlockObjectiveValue value = evaluate(state);
    SpectralIncrement constant_gradient;
    switch (block_) {
        case LocalSldBlock::full:
            constant_gradient =
                full_objective.signed_constant_ratio_gradient(state);
            break;
        case LocalSldBlock::selected_closed:
            constant_gradient =
                selected_objective.signed_constant_ratio_gradient(state);
            break;
        case LocalSldBlock::complement_closed:
            constant_gradient =
                complement_objective.signed_constant_ratio_gradient(state);
            break;
        case LocalSldBlock::mixed:
            constant_gradient =
                full_objective.signed_constant_ratio_gradient(state);
            add_scaled(
                constant_gradient,
                selected_objective.signed_constant_ratio_gradient(state),
                -1.0L);
            add_scaled(
                constant_gradient,
                complement_objective.signed_constant_ratio_gradient(state),
                -1.0L);
            break;
    }
    if (!(full.energy > 0.0L) || !(full.enstrophy > 0.0L) ||
        !(full.palinstrophy > 0.0L)) {
        return SpectralIncrement(state.waves.size());
    }
    const SpectralReal x = value.normalized_stretching;
    const SpectralReal x2 = x * x;
    const SpectralReal x4 = x2 * x2;
    const SpectralReal shape_derivative =
        4.0L * x2 * (3.0L - x4) /
        ((1.0L + x4) * (1.0L + x4));
    const SpectralReal stretching_scale =
        std::pow(full.energy, 0.25L) *
        std::pow(full.enstrophy, 0.25L) *
        full.palinstrophy;
    SpectralIncrement x_gradient =
        full_objective.signed_stretching_gradient(state);
    scale(x_gradient, 1.0L / stretching_scale);
    add_scaled(
        x_gradient, state.velocity,
        -0.5L * x / full.energy);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    add_scaled(
        x_gradient, au,
        -0.5L * x / full.enstrophy);
    add_scaled(
        x_gradient, laplacian_weight(state, au),
        -2.0L * x / full.palinstrophy);

    scale(constant_gradient, value.common_shape_factor);
    add_scaled(
        constant_gradient, x_gradient,
        value.block_constant_ratio * shape_derivative);
    return constant_gradient;
}

}  // namespace lemma
