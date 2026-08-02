#include "local_sld_projective_open_power_objective.hpp"

#include "projective_advection_decomposition.hpp"
#include "projective_quartic_diagonal_kernel.hpp"

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
            "projective open-power gradient layout mismatch");
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
            "projective open-power Laplacian layout mismatch");
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

LocalSldProjectiveOpenPowerObjective::
LocalSldProjectiveOpenPowerObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads)
    : dynamics_(dynamics),
      selection_(selection),
      core_maximum_height_(core_maximum_height),
      threads_(threads),
      core_(ProjectiveCoreFamily::through_maximum_height(
          core_maximum_height)) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective open-power threads must be 1..256");
    }
}

LocalSldProjectiveOpenPowerObjectiveValue
LocalSldProjectiveOpenPowerObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    LocalSldProjectiveOpenPowerObjectiveValue result;
    result.core_maximum_height = core_maximum_height_;
    result.fixed_core_shape_count = core_.size();
    result.enstrophy = selected.enstrophy;
    result.palinstrophy = selected.palinstrophy;
    result.full_stretching = full.signed_stretching;
    result.selected_bracket = selected.signed_two_entry_bracket;
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    const auto& aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, core_);
    result.fixed_core_bracket =
        ProjectiveQuarticDiagonalKernel::evaluate(
            state, aggregate, result.enstrophy,
            result.palinstrophy, false, threads_).bracket;
    result.open_bracket =
        result.selected_bracket - result.fixed_core_bracket;
    result.signed_open_power_one =
        result.open_bracket * result.full_stretching /
        (result.enstrophy * result.enstrophy *
         result.palinstrophy * result.palinstrophy);
    result.absolute_open_power_one =
        std::abs(result.signed_open_power_one);
    result.squared_open_power_one =
        result.signed_open_power_one *
        result.signed_open_power_one;
    result.finite = std::isfinite(result.squared_open_power_one);
    return result;
}

SpectralIncrement LocalSldProjectiveOpenPowerObjective::gradient(
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
    const auto& aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, core_);
    const ProjectiveQuarticDiagonalMoment core_moment =
        ProjectiveQuarticDiagonalKernel::evaluate(
            state, aggregate, selected.enstrophy,
            selected.palinstrophy, true, threads_);
    const SpectralReal open =
        selected.signed_two_entry_bracket - core_moment.bracket;
    const SpectralReal denominator =
        selected.enstrophy * selected.enstrophy *
        selected.palinstrophy * selected.palinstrophy;
    const SpectralReal signed_power =
        open * full.signed_stretching / denominator;

    SpectralIncrement open_gradient =
        selected_objective.two_entry_bracket_gradient(state);
    add_scaled(open_gradient, core_moment.gradient, -1.0L);
    result = open_gradient;
    scale(result, full.signed_stretching / denominator);
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        open / denominator);
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
