#include "local_sld_projective_height_outer_power_objective.hpp"

#include "projective_advection_decomposition.hpp"
#include "projective_height_shell_partition.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
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

struct OuterGraph {
    std::size_t shell_count = 0;
    std::size_t active_shell_count = 0;
    SpectralReal h1_sum = 0.0L;
    SpectralIncrement gradient;
};

OuterGraph build_outer_graph(
    const SpectralState& state,
    TriadSelection selection,
    int threads,
    bool compute_gradient) {
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    const auto shells = ProjectiveHeightShellPartition::build(groups);
    OuterGraph result;
    result.shell_count = shells.size();
    if (compute_gradient) {
        result.gradient.resize(state.waves.size());
    }
    for (const auto& shell : shells) {
        if (shell.group_indices.empty()) {
            continue;
        }
        ++result.active_shell_count;
        const SpectralIncrement b =
            ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
                state, groups, shell.group_indices,
                state.velocity, state.velocity, threads);
        SpectralIncrement ab = laplacian_weight(state, b);
        result.h1_sum += pairing(b, ab);
        if (!compute_gradient) {
            continue;
        }
        scale(ab, 2.0L);
        const auto& aggregate =
            ProjectiveAdvectionDecomposition::aggregate_family(
                state, selection, shell.primitive_shapes);
        add_scaled(
            result.gradient,
            ProjectiveAdvectionDecomposition::vjp(
                state, aggregate.front(), ab, threads),
            1.0L);
    }
    return result;
}

}  // namespace

LocalSldProjectiveHeightOuterPowerObjective::
LocalSldProjectiveHeightOuterPowerObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    int threads)
    : dynamics_(dynamics), selection_(selection), threads_(threads) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height outer-power threads must be 1..256");
    }
}

LocalSldProjectiveHeightOuterPowerObjectiveValue
LocalSldProjectiveHeightOuterPowerObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    const OuterGraph outer = build_outer_graph(
        state, selection_, threads_, false);
    LocalSldProjectiveHeightOuterPowerObjectiveValue result;
    result.height_shell_count = outer.shell_count;
    result.active_height_shell_count = outer.active_shell_count;
    result.enstrophy = selected.enstrophy;
    result.palinstrophy = selected.palinstrophy;
    result.full_stretching = full.signed_stretching;
    result.diagonal_outer_h1_sum = outer.h1_sum;
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal denominator =
        result.enstrophy * result.enstrophy *
        result.palinstrophy * result.palinstrophy;
    result.signed_outer_power_one =
        result.diagonal_outer_h1_sum * result.full_stretching /
        denominator;
    result.absolute_outer_power_one =
        std::abs(result.signed_outer_power_one);
    result.squared_outer_power_one =
        result.signed_outer_power_one *
        result.signed_outer_power_one;
    result.finite = std::isfinite(result.squared_outer_power_one);
    return result;
}

SpectralIncrement
LocalSldProjectiveHeightOuterPowerObjective::gradient(
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
    const OuterGraph outer = build_outer_graph(
        state, selection_, threads_, true);
    const SpectralReal denominator =
        selected.enstrophy * selected.enstrophy *
        selected.palinstrophy * selected.palinstrophy;
    const SpectralReal signed_power =
        outer.h1_sum * full.signed_stretching / denominator;
    result = outer.gradient;
    scale(result, full.signed_stretching / denominator);
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        outer.h1_sum / denominator);
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
