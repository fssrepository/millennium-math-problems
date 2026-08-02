#include "local_sld_projective_height_power_objective.hpp"

#include "projective_advection_decomposition.hpp"
#include "projective_quartic_diagonal_kernel.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

std::vector<Shape> shell_shapes(
    const std::vector<ProjectiveInteractionGroup>& groups,
    SpectralInteger lower,
    SpectralInteger upper) {
    std::vector<Shape> result;
    for (const auto& group : groups) {
        const SpectralInteger height =
            group.primitive_squared_lengths[2];
        if (height > lower && height <= upper) {
            result.push_back(group.primitive_squared_lengths);
        }
    }
    return result;
}

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

LocalSldProjectiveHeightPowerObjective::
LocalSldProjectiveHeightPowerObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads)
    : dynamics_(dynamics),
      selection_(selection),
      core_maximum_height_(core_maximum_height),
      shell_maximum_height_(2 * core_maximum_height),
      threads_(threads) {
    if (core_maximum_height < 1) {
        throw std::invalid_argument(
            "projective height-power core height must be positive");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height-power threads must be 1..256");
    }
}

LocalSldProjectiveHeightPowerObjectiveValue
LocalSldProjectiveHeightPowerObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    const std::vector<Shape> shapes = shell_shapes(
        groups, core_maximum_height_, shell_maximum_height_);
    LocalSldProjectiveHeightPowerObjectiveValue result;
    result.core_maximum_height = core_maximum_height_;
    result.shell_maximum_height = shell_maximum_height_;
    result.shell_shape_count = shapes.size();
    result.enstrophy = selected.enstrophy;
    result.palinstrophy = selected.palinstrophy;
    result.full_stretching = full.signed_stretching;
    if (shapes.empty() || !(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    const auto& aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, shapes);
    result.shell_internal_bracket =
        ProjectiveQuarticDiagonalKernel::evaluate(
            state, aggregate, result.enstrophy,
            result.palinstrophy, false, threads_).bracket;
    const SpectralReal denominator =
        result.enstrophy * result.enstrophy *
        result.palinstrophy * result.palinstrophy;
    result.signed_shell_power_one =
        result.shell_internal_bracket * result.full_stretching /
        denominator;
    result.absolute_shell_power_one =
        std::abs(result.signed_shell_power_one);
    result.squared_shell_power_one =
        result.signed_shell_power_one *
        result.signed_shell_power_one;
    result.finite = std::isfinite(result.squared_shell_power_one);
    return result;
}

SpectralIncrement LocalSldProjectiveHeightPowerObjective::gradient(
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
    const std::vector<Shape> shapes = shell_shapes(
        groups, core_maximum_height_, shell_maximum_height_);
    if (shapes.empty()) {
        return result;
    }
    const auto& aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, shapes);
    const ProjectiveQuarticDiagonalMoment shell_moment =
        ProjectiveQuarticDiagonalKernel::evaluate(
            state, aggregate, selected.enstrophy,
            selected.palinstrophy, true, threads_);
    const SpectralReal denominator =
        selected.enstrophy * selected.enstrophy *
        selected.palinstrophy * selected.palinstrophy;
    const SpectralReal signed_power =
        shell_moment.bracket * full.signed_stretching /
        denominator;
    result = shell_moment.gradient;
    scale(result, full.signed_stretching / denominator);
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        shell_moment.bracket / denominator);
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
