#include "local_sld_projective_cross_power_objective.hpp"
#include "projective_quartic_diagonal_kernel.hpp"

#include <algorithm>
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

void add_scaled(
    SpectralIncrement& target,
    const SpectralIncrement& source,
    SpectralReal factor) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "projective cross-power gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] +=
                factor * source[mode][coordinate];
        }
    }
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective cross-power Laplacian layout mismatch");
    }
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

}  // namespace

LocalSldProjectiveCrossPowerObjective::
LocalSldProjectiveCrossPowerObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    int threads)
    : dynamics_(dynamics), selection_(selection), threads_(threads) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective cross-power threads must be 1..256");
    }
}

SpectralReal LocalSldProjectiveCrossPowerObjective::diagonal_bracket(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L)) {
        return 0.0L;
    }
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    return ProjectiveQuarticDiagonalKernel::evaluate(
        state, groups, selected.enstrophy,
        selected.palinstrophy, false, threads_).bracket;
}

SpectralIncrement
LocalSldProjectiveCrossPowerObjective::diagonal_bracket_gradient(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L)) {
        return SpectralIncrement(state.waves.size());
    }
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    return ProjectiveQuarticDiagonalKernel::evaluate(
        state, groups, selected.enstrophy,
        selected.palinstrophy, true, threads_).gradient;
}

LocalSldProjectiveCrossPowerObjectiveValue
LocalSldProjectiveCrossPowerObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    LocalSldProjectiveCrossPowerObjectiveValue result;
    result.projective_shape_count =
        ProjectiveAdvectionDecomposition::group(
            state, selection_).size();
    result.enstrophy = selected.enstrophy;
    result.palinstrophy = selected.palinstrophy;
    result.full_stretching = full.signed_stretching;
    result.full_bracket = selected.signed_two_entry_bracket;
    result.diagonal_bracket = diagonal_bracket(state);
    result.cross_bracket =
        result.full_bracket - result.diagonal_bracket;
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    result.signed_cross_power_one =
        result.cross_bracket * result.full_stretching /
        (result.enstrophy * result.enstrophy *
         result.palinstrophy * result.palinstrophy);
    result.absolute_cross_power_one =
        std::abs(result.signed_cross_power_one);
    result.squared_cross_power_one =
        result.signed_cross_power_one *
        result.signed_cross_power_one;
    result.finite = std::isfinite(result.squared_cross_power_one);
    return result;
}

SpectralIncrement LocalSldProjectiveCrossPowerObjective::gradient(
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
        !(selected.palinstrophy > 0.0L)) {
        return result;
    }
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    const ProjectiveQuarticDiagonalMoment diagonal_moment =
        ProjectiveQuarticDiagonalKernel::evaluate(
            state, groups, selected.enstrophy,
            selected.palinstrophy, true, threads_);
    const SpectralReal diagonal = diagonal_moment.bracket;
    const SpectralReal cross =
        selected.signed_two_entry_bracket - diagonal;
    const SpectralReal denominator =
        selected.enstrophy * selected.enstrophy *
        selected.palinstrophy * selected.palinstrophy;
    const SpectralReal signed_power =
        cross * full.signed_stretching / denominator;

    SpectralIncrement cross_gradient =
        selected_objective.two_entry_bracket_gradient(state);
    add_scaled(cross_gradient, diagonal_moment.gradient, -1.0L);
    result = cross_gradient;
    scale(
        result,
        full.signed_stretching / denominator);
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        cross / denominator);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    add_scaled(
        result, au,
        -4.0L * signed_power / selected.enstrophy);
    add_scaled(
        result, laplacian_weight(state, au),
        -4.0L * signed_power / selected.palinstrophy);
    scale(result, 2.0L * signed_power);
    return result;
}

}  // namespace lemma
