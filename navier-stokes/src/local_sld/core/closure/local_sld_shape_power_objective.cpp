#include "local_sld_shape_power_objective.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

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
        throw std::invalid_argument(
            "shape-power gradient layout mismatch");
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

SpectralReal integer_power(SpectralReal value, int power) {
    SpectralReal result = 1.0L;
    for (int index = 0; index < power; ++index) {
        result *= value;
    }
    return result;
}

}  // namespace

LocalSldShapePowerObjective::LocalSldShapePowerObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    int power)
    : dynamics_(dynamics), selection_(selection), power_(power) {
    if (power < 0 || power > 3) {
        throw std::invalid_argument(
            "shape-power objective supports integer powers 0 through 3");
    }
}

LocalSldShapePowerObjectiveValue
LocalSldShapePowerObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    LocalSldShapePowerObjectiveValue result;
    result.power = power_;
    result.bracket_constant_ratio = selected.signed_constant_ratio;
    result.normalized_stretching = full.normalized_stretching_ratio;
    const SpectralReal signed_product =
        result.bracket_constant_ratio *
        integer_power(result.normalized_stretching, power_);
    result.absolute_power_product = std::abs(signed_product);
    result.squared_power_product = signed_product * signed_product;
    result.finite = selected.finite && full.finite &&
        std::isfinite(result.squared_power_product);
    return result;
}

SpectralIncrement LocalSldShapePowerObjective::gradient(
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
    if (!(full.energy > 0.0L) || !(full.enstrophy > 0.0L) ||
        !(full.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal c = selected.signed_constant_ratio;
    const SpectralReal x = full.normalized_stretching_ratio;
    const SpectralReal x_power = integer_power(x, power_);
    const SpectralReal x_power2 = x_power * x_power;
    result = selected_objective.signed_constant_ratio_gradient(state);
    scale(result, 2.0L * c * x_power2);
    if (power_ == 0) {
        return result;
    }

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
    const SpectralReal x_derivative_power = integer_power(
        x, 2 * power_ - 1);
    add_scaled(
        result, x_gradient,
        2.0L * static_cast<SpectralReal>(power_) * c * c *
            x_derivative_power);
    return result;
}

}  // namespace lemma
