#include "local_quartic_closure_adversary.hpp"

#include "initial_sobolev_constraint.hpp"
#include "local_sld_trajectory_adjoint.hpp"
#include "parallel_executor.hpp"
#include "spectral_adjoint.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_objective.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>

namespace lemma {
namespace {

TriadSelection closure_selection(const std::string& name) {
    if (name == "local") {
        return TriadPartition::local;
    }
    if (name == "doubling-family") {
        return TriadSelection::local_equal_low_doubling();
    }
    if (name == "doubling-remainder") {
        return TriadSelection::local_without_equal_low_doubling();
    }
    throw std::invalid_argument(
        "closure selection must be local, doubling-family, or doubling-remainder");
}

bool is_common_block_objective(const std::string& objective) {
    return objective == "block-ratio" || objective == "mixed-ratio";
}

bool is_frozen_trajectory_objective(const std::string& objective) {
    return objective == "terminal-sld-ratio" ||
        objective == "maximum-sld-ratio";
}

LocalSldBlock block_for_objective(const std::string& objective) {
    return objective == "mixed-ratio"
        ? LocalSldBlock::mixed
        : LocalSldBlock::selected_closed;
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

SpectralReal state_distance_squared(const SpectralState& left,
                                    const SpectralState& right) {
    if (left.waves != right.waves ||
        left.velocity.size() != right.velocity.size()) {
        return std::numeric_limits<SpectralReal>::infinity();
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            result += std::norm(
                left.velocity[mode][component] -
                right.velocity[mode][component]);
        }
    }
    return result;
}

SpectralState make_start(
    int cutoff, int restart, std::uint64_t seed,
    const SpectralState* previous_winner,
    const InitialSobolevConstraint& sobolev) {
    std::mt19937_64 generator(seed);
    SpectralState state;
    if (previous_winner != nullptr && restart == 0) {
        state = SpectralStateFactory::lift(
            *previous_winner, cutoff, generator);
    } else if (previous_winner != nullptr && restart % 3 == 0) {
        state = SpectralStateFactory::lift(
            *previous_winner, cutoff, generator);
        state = SpectralStateFactory::mutate(
            state, 0.08L + 0.02L * static_cast<SpectralReal>(restart % 5),
            generator, restart % 2 == 0);
    } else if (restart % 2 == 0) {
        const SpectralReal decay = 0.55L +
            0.08L * static_cast<SpectralReal>(restart % 5);
        state = SpectralStateFactory::analytic(cutoff, seed, decay);
    } else {
        state = SpectralStateFactory::random(cutoff, generator);
    }
    sobolev.retract(state, 1.0L);
    return state;
}

SpectralReal fitted_slope(
    const std::vector<LocalQuarticClosureCutoffResult>& rows) {
    SpectralReal n = 0.0L;
    SpectralReal sx = 0.0L;
    SpectralReal sy = 0.0L;
    SpectralReal sxx = 0.0L;
    SpectralReal sxy = 0.0L;
    for (const auto& row : rows) {
        if (row.cutoff <= 1 ||
            !(row.winner.objective > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
        const SpectralReal y = std::log(
            row.winner.objective);
        n += 1.0L;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    const SpectralReal denominator = n * sxx - sx * sx;
    return n >= 2.0L && std::abs(denominator) > 1e-30L
        ? (n * sxy - sx * sy) / denominator
        : 0.0L;
}

}  // namespace

LocalQuarticClosureRestartResult
LocalQuarticClosureAdversary::maximize(
    const SpectralState& initial,
    const LocalQuarticClosureAdversaryOptions& options,
    int restart, std::uint64_t seed, bool warm_continuation) {
    SpectralGalerkin galerkin;
    galerkin.configure(options.backend, 1);
    const SpectralDynamics dynamics(galerkin);
    const SpectralObjective spectral_objective(dynamics);
    const SpectralAdjoint adjoint(dynamics, spectral_objective);
    const GradientAdversary adversary(
        dynamics, spectral_objective, adjoint);
    GradientSearchOptions search;
    search.iterations = options.iterations;
    search.line_search_steps = options.line_search_steps;
    search.trajectory_steps = is_frozen_trajectory_objective(
        options.objective)
        ? options.trajectory_steps
        : 0;
    search.viscosity = options.viscosity;
    search.time_step = options.time_step;
    search.initial_step = options.initial_step;
    search.objective = options.objective == "closure-ratio"
        ? "local-closure-ratio"
        : (options.objective == "signed-closure-ratio"
               ? "local-signed-closure-ratio"
               : (options.objective == "block-ratio"
                      ? "local-sld-block-ratio"
                      : (options.objective == "mixed-ratio"
                             ? "local-sld-mixed-ratio"
                             : (options.objective == "terminal-sld-ratio"
                                    ? "local-frozen-terminal-sld-ratio"
                                    : (options.objective ==
                                               "maximum-sld-ratio"
                                           ? "local-frozen-maximum-sld-ratio"
                                           : "local-sld-ratio")))));
    search.method = options.method;
    search.lbfgs_history = options.lbfgs_history;
    search.sobolev_order = options.sobolev_order;
    search.sobolev_cap = options.sobolev_cap;
    search.closure_selection = closure_selection(options.selection);
    const GradientSearchResult optimized = adversary.maximize_q(
        initial, search);
    const LocalQuarticClosureObjective closure(
        dynamics, search.closure_selection);
    const LocalQuarticClosureObjectiveValue initial_value =
        closure.evaluate(initial);
    LocalQuarticClosureRestartResult result;
    result.state = optimized.state;
    result.value = closure.evaluate(result.state);
    result.common_block_objective =
        is_common_block_objective(options.objective);
    if (result.common_block_objective) {
        result.common_block_value = LocalSldBlockObjective(
            dynamics, search.closure_selection,
            block_for_objective(options.objective))
            .evaluate(result.state);
    }
    if (is_frozen_trajectory_objective(options.objective)) {
        const LocalSldTrajectoryAdjoint trajectory(
            dynamics, search.closure_selection);
        const LocalSldTrajectoryValue coarse =
            options.objective == "maximum-sld-ratio"
            ? trajectory.maximum_value(
                  result.state, options.viscosity, options.time_step,
                  options.trajectory_steps)
            : trajectory.terminal_value(
                  result.state, options.viscosity, options.time_step,
                  options.trajectory_steps);
        const LocalSldTrajectoryValue refined =
            options.objective == "maximum-sld-ratio"
            ? trajectory.maximum_value(
                  result.state, options.viscosity,
                  options.time_step / 2.0L,
                  options.trajectory_steps * 2)
            : trajectory.terminal_value(
                  result.state, options.viscosity,
                  options.time_step / 2.0L,
                  options.trajectory_steps * 2);
        result.refined_objective = refined.terminal_ratio;
        result.time_step_relative_error = std::abs(
            refined.terminal_ratio - coarse.terminal_ratio) /
            std::max({std::abs(refined.terminal_ratio),
                      std::abs(coarse.terminal_ratio), 1e-30L});
        result.frozen_initial_frequency = coarse.initial_frequency;
        result.frozen_initial_ep_shift = coarse.initial_ep_shift;
        result.objective_step = coarse.steps;
        result.refined_objective_step = refined.steps;
    }
    result.initial_objective = optimized.initial_objective;
    result.objective = optimized.objective;
    if (!is_frozen_trajectory_objective(options.objective)) {
        result.objective_step = optimized.objective_step;
    }
    result.initial_constant_ratio = initial_value.constant_ratio;
    result.final_projected_gradient_norm =
        optimized.final_projected_gradient_norm;
    result.sobolev_value = optimized.final_sobolev_value;
    result.seed = seed;
    result.restart = restart;
    result.accepted_steps = optimized.accepted_steps;
    result.evaluations = optimized.trajectory_evaluations;
    result.warm_continuation = warm_continuation;
    return result;
}

LocalQuarticClosureAdversaryReport LocalQuarticClosureEnsemble::scan(
    const LocalQuarticClosureAdversaryOptions& options) {
    if (options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 8 || options.restarts < 1 ||
        options.restarts > 1000 || options.workers < 1 ||
        options.workers > 256 || options.iterations < 1 ||
        options.line_search_steps < 1 ||
        !(options.initial_step > 0.0L) ||
        (options.method != "steepest" && options.method != "lbfgs") ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft") ||
        (options.objective != "closure-ratio" &&
         options.objective != "signed-closure-ratio" &&
         options.objective != "sld-ratio" &&
         options.objective != "block-ratio" &&
         options.objective != "mixed-ratio" &&
         options.objective != "terminal-sld-ratio" &&
         options.objective != "maximum-sld-ratio") ||
        (options.selection != "local" &&
         options.selection != "doubling-family" &&
         options.selection != "doubling-remainder") ||
        (is_common_block_objective(options.objective) &&
         options.selection == "local") ||
        (is_frozen_trajectory_objective(options.objective) &&
         (options.trajectory_steps < 1 || !(options.viscosity > 0.0L) ||
          !(options.time_step > 0.0L)))) {
        throw std::invalid_argument(
            "invalid local quartic closure adversary options");
    }
    const InitialSobolevConstraint sobolev(
        options.sobolev_order, options.sobolev_cap);
    const ParallelExecutor executor(options.workers);
    LocalQuarticClosureAdversaryReport report;
    report.workers = executor.threads();
    report.restarts = options.restarts;
    report.iterations = options.iterations;
    report.trajectory_steps = is_frozen_trajectory_objective(
        options.objective)
        ? options.trajectory_steps
        : 0;
    report.objective = options.objective;
    report.backend = options.backend;
    report.viscosity = options.viscosity;
    report.time_step = options.time_step;
    report.sobolev_order = options.sobolev_order;
    report.sobolev_cap = options.sobolev_cap;
    SpectralState previous_winner;
    bool has_previous_winner = false;
    if (!options.warm_state_path.empty()) {
        previous_winner = SpectralStateReader::read_tsv(
            options.warm_state_path);
        const int warm_cutoff = SpectralStateOps::cutoff(previous_winner);
        if (warm_cutoff != options.minimum_cutoff - 1 &&
            warm_cutoff != options.minimum_cutoff) {
            throw std::invalid_argument(
                "closure warm state cutoff must equal min-cutoff or min-cutoff minus one");
        }
        sobolev.retract(previous_winner, 1.0L);
        has_previous_winner = true;
    }

    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        std::vector<SpectralState> starts(
            static_cast<std::size_t>(options.restarts));
        std::vector<std::uint64_t> seeds(
            static_cast<std::size_t>(options.restarts));
        for (int restart = 0; restart < options.restarts; ++restart) {
            const std::uint64_t seed = splitmix64(
                options.seed ^
                splitmix64(static_cast<std::uint64_t>(cutoff)) ^
                splitmix64(static_cast<std::uint64_t>(restart + 1)));
            seeds[static_cast<std::size_t>(restart)] = seed;
            starts[static_cast<std::size_t>(restart)] = make_start(
                cutoff, restart, seed,
                has_previous_winner ? &previous_winner : nullptr,
                sobolev);
        }

        LocalQuarticClosureCutoffResult row;
        row.cutoff = cutoff;
        if (has_previous_winner) {
            SpectralGalerkin galerkin;
            galerkin.configure(options.backend, 1);
            const SpectralDynamics dynamics(galerkin);
            const LocalQuarticClosureObjectiveValue warm_value =
                LocalQuarticClosureObjective(
                    dynamics, closure_selection(options.selection))
                    .evaluate(starts.front());
            row.warm_lift_constant_ratio = warm_value.constant_ratio;
            row.warm_lift_objective = options.objective == "closure-ratio"
                ? warm_value.squared_constant_ratio
                : (options.objective == "signed-closure-ratio"
                       ? warm_value.signed_constant_ratio
                       : (options.objective == "terminal-sld-ratio"
                              ? LocalSldTrajectoryAdjoint(
                                    dynamics,
                                    closure_selection(options.selection))
                                    .terminal_value(
                                        starts.front(), options.viscosity,
                                        options.time_step,
                                        options.trajectory_steps)
                                    .terminal_ratio
                       : (options.objective == "maximum-sld-ratio"
                              ? LocalSldTrajectoryAdjoint(
                                    dynamics,
                                    closure_selection(options.selection))
                                    .maximum_value(
                                        starts.front(), options.viscosity,
                                        options.time_step,
                                        options.trajectory_steps)
                                    .terminal_ratio
                       : (is_common_block_objective(options.objective)
                              ? LocalSldBlockObjective(
                                    dynamics,
                                    closure_selection(options.selection),
                                    block_for_objective(options.objective))
                                    .evaluate(starts.front())
                                    .block_sld_ratio
                              : warm_value.signed_local_sld_ratio))));
        }
        std::vector<LocalQuarticClosureRestartResult> results(
            static_cast<std::size_t>(options.restarts));
        executor.for_each(results.size(), [&](std::size_t restart) {
            results[restart] = LocalQuarticClosureAdversary::maximize(
                starts[restart], options, static_cast<int>(restart),
                seeds[restart], has_previous_winner && restart == 0);
        });
        row.restart_constant_ratios.reserve(results.size());
        for (auto& candidate : results) {
            row.restart_constant_ratios.push_back(
                candidate.value.constant_ratio);
            row.restart_objectives.push_back(candidate.objective);
            if (!row.winner.value.finite ||
                candidate.objective > row.winner.objective) {
                row.winner = std::move(candidate);
            }
        }
        row.improvement_factor = row.winner.initial_constant_ratio > 0.0L
            ? row.winner.value.constant_ratio /
                  row.winner.initial_constant_ratio
            : 0.0L;
        row.objective_gain = row.winner.objective -
            row.winner.initial_objective;
        if (has_previous_winner) {
            if (SpectralStateOps::cutoff(previous_winner) == cutoff) {
                row.projection_residual = std::sqrt(
                    state_distance_squared(
                        row.winner.state, previous_winner));
            } else {
                const SpectralState projected =
                    SpectralStateFactory::project(
                        row.winner.state, cutoff - 1);
                row.projection_residual = std::sqrt(
                    state_distance_squared(projected, previous_winner));
            }
        }
        report.maximum_constant_ratio = std::max(
            report.maximum_constant_ratio,
            row.winner.value.constant_ratio);
        report.maximum_objective = std::max(
            report.maximum_objective, row.winner.objective);
        previous_winner = row.winner.state;
        has_previous_winner = true;
        report.rows.push_back(std::move(row));
    }
    report.fitted_cutoff_slope = fitted_slope(report.rows);
    return report;
}

}  // namespace lemma
