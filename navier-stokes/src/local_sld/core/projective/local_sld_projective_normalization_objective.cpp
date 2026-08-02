#include "local_sld_projective_normalization_objective.hpp"

#include <cmath>

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
    TriadSelection selection)
    : dynamics_(dynamics), selection_(selection) {}

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
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal bracket =
        1.5L * result.selected_stretching *
        result.selected_palinstrophy_cross / result.palinstrophy;
    result.palinstrophy_normalization_power_one =
        std::abs(result.full_stretching * bracket) /
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
        selected.signed_stretching == 0.0L ||
        selected.palinstrophy_cross == 0.0L) {
        return result;
    }
    const SpectralReal value = evaluate(state)
        .squared_palinstrophy_normalization_power_one;
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        2.0L * value / full.signed_stretching);
    add_scaled(
        result,
        selected_objective.signed_stretching_gradient(state),
        2.0L * value / selected.signed_stretching);
    add_scaled(
        result,
        selected_objective.palinstrophy_cross_gradient(state),
        2.0L * value / selected.palinstrophy_cross);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    add_scaled(result, au, -8.0L * value / selected.enstrophy);
    add_scaled(
        result, laplacian_weight(state, au),
        -12.0L * value / selected.palinstrophy);
    return result;
}

}  // namespace lemma
