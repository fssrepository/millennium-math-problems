#include "local_sld_projective_normalization_objective.hpp"

#include "projective_advection_decomposition.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

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

void add_scaled(
    SpectralIncrement& target,
    const SpectralIncrement& source,
    SpectralReal factor) {
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] +=
                factor * source[mode][coordinate];
        }
    }
}

}  // namespace

LocalSldProjectiveNormalizationObjective::
LocalSldProjectiveNormalizationObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads)
    : dynamics_(dynamics),
      selection_(selection),
      core_maximum_height_(core_maximum_height),
      threads_(threads),
      core_(core_maximum_height > 0
          ? ProjectiveCoreFamily::through_maximum_height(
                core_maximum_height)
          : std::vector<ProjectivePrimitiveSignature>{}) {
    if (threads < 1 || threads > 256 || core_maximum_height < 0) {
        throw std::invalid_argument(
            "projective normalization requires nonnegative core height and threads 1..256");
    }
}

LocalSldProjectiveNormalizationObjectiveValue
LocalSldProjectiveNormalizationObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    LocalSldProjectiveNormalizationObjectiveValue result;
    result.enstrophy = selected.enstrophy;
    result.palinstrophy = selected.palinstrophy;
    result.selected_stretching = selected.signed_stretching;
    result.selected_palinstrophy_cross = selected.palinstrophy_cross;
    result.full_stretching = full.signed_stretching;
    result.core_maximum_height = core_maximum_height_;
    result.fixed_core_shape_count = core_.size();
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    if (!core_.empty()) {
        const auto& aggregate =
            ProjectiveAdvectionDecomposition::aggregate_family(
                state, selection_, core_);
        const SpectralIncrement core_b =
            ProjectiveAdvectionDecomposition::evaluate_bilinear(
                state, aggregate.front(),
                state.velocity, state.velocity, threads_);
        const SpectralIncrement au = laplacian_weight(
            state, state.velocity);
        const SpectralIncrement core_ab = laplacian_weight(state, core_b);
        SpectralReal core_stretching = 0.0L;
        SpectralReal core_cross = 0.0L;
        for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
            core_stretching += std::real(dot_hermitian(
                au[mode], core_b[mode]));
            core_cross += std::real(dot_hermitian(
                core_ab[mode], au[mode]));
        }
        result.fixed_core_stretching = core_stretching;
        result.fixed_core_palinstrophy_cross = core_cross;
    }
    result.open_palinstrophy_normalization =
        1.5L *
        (result.selected_stretching *
             result.selected_palinstrophy_cross -
         result.fixed_core_stretching *
             result.fixed_core_palinstrophy_cross) /
        result.palinstrophy;
    result.palinstrophy_normalization_power_one =
        std::abs(
            result.full_stretching *
            result.open_palinstrophy_normalization) /
        (result.enstrophy * result.enstrophy *
         result.palinstrophy * result.palinstrophy);
    result.squared_palinstrophy_normalization_power_one =
        result.palinstrophy_normalization_power_one *
        result.palinstrophy_normalization_power_one;
    result.finite = std::isfinite(
        result.squared_palinstrophy_normalization_power_one);
    return result;
}

SpectralIncrement LocalSldProjectiveNormalizationObjective::gradient(
    const SpectralState& state) const {
    const LocalQuarticClosureObjective selected_objective(
        dynamics_, selection_);
    const LocalQuarticClosureObjective full_objective(
        dynamics_, TriadPartition::local);
    const LocalQuarticClosureObjectiveValue selected =
        selected_objective.evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        full_objective.evaluate(state);
    SpectralIncrement result(state.waves.size());
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L) ||
        full.signed_stretching == 0.0L ||
        selected.signed_stretching == 0.0L) {
        return result;
    }
    const LocalSldProjectiveNormalizationObjectiveValue objective_value =
        evaluate(state);
    const SpectralReal value = objective_value
        .squared_palinstrophy_normalization_power_one;
    const SpectralReal open_product =
        selected.signed_stretching * selected.palinstrophy_cross -
        objective_value.fixed_core_stretching *
            objective_value.fixed_core_palinstrophy_cross;
    if (open_product == 0.0L) {
        return result;
    }
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        2.0L * value / full.signed_stretching);
    SpectralIncrement open_product_gradient =
        selected_objective.signed_stretching_gradient(state);
    for (ComplexVector& mode : open_product_gradient) {
        for (SpectralComplex& component : mode) {
            component *= selected.palinstrophy_cross;
        }
    }
    add_scaled(
        open_product_gradient,
        selected_objective.palinstrophy_cross_gradient(state),
        selected.signed_stretching);
    if (!core_.empty()) {
        const auto& aggregate =
            ProjectiveAdvectionDecomposition::aggregate_family(
                state, selection_, core_);
        const auto& core_group = aggregate.front();
        const SpectralIncrement core_b =
            ProjectiveAdvectionDecomposition::evaluate_bilinear(
                state, core_group,
                state.velocity, state.velocity, threads_);
        const SpectralIncrement au = laplacian_weight(
            state, state.velocity);
        const SpectralIncrement core_ab = laplacian_weight(state, core_b);
        SpectralIncrement core_stretching_gradient = core_ab;
        add_scaled(
            core_stretching_gradient,
            ProjectiveAdvectionDecomposition::vjp(
                state, core_group, au, threads_),
            1.0L);
        SpectralIncrement core_cross_gradient =
            laplacian_weight(state, core_ab);
        add_scaled(
            core_cross_gradient,
            ProjectiveAdvectionDecomposition::vjp(
                state, core_group,
                laplacian_weight(state, au), threads_),
            1.0L);
        add_scaled(
            open_product_gradient,
            core_stretching_gradient,
            -objective_value.fixed_core_palinstrophy_cross);
        add_scaled(
            open_product_gradient,
            core_cross_gradient,
            -objective_value.fixed_core_stretching);
    }
    add_scaled(
        result, open_product_gradient,
        2.0L * value / open_product);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    add_scaled(result, au, -8.0L * value / selected.enstrophy);
    add_scaled(
        result, laplacian_weight(state, au),
        -12.0L * value / selected.palinstrophy);
    return result;
}

}  // namespace lemma
