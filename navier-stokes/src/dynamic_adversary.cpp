#include "dynamic_adversary.hpp"

#include "initial_sobolev_constraint.hpp"
#include "parallel_executor.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>

namespace lemma {
namespace {

bool collects_partition(const std::string& objective) {
    return objective == "critical-local-integral" ||
           objective == "critical-local-increase" ||
           objective == "critical-nonlocal-integral" ||
           objective == "critical-near-nonlocal-integral" ||
           objective == "critical-far-nonlocal-integral" ||
           objective == "critical-gap-tail-integral";
}

SpectralState move_to_cutoff(const SpectralState& source, int target_cutoff,
                             std::mt19937_64& generator) {
    const int source_cutoff = SpectralStateOps::cutoff(source);
    return source_cutoff <= target_cutoff
               ? SpectralStateFactory::lift(source, target_cutoff, generator)
               : SpectralStateFactory::project(source, target_cutoff);
}

}  // namespace

DynamicAdversary::DynamicAdversary(std::string backend, int compute_threads)
    : dynamics_(galerkin_),
      objective_(dynamics_),
      trajectory_(galerkin_, dynamics_, objective_),
      adjoint_(dynamics_, objective_),
      gradient_(dynamics_, objective_, adjoint_) {
    galerkin_.configure(backend, compute_threads);
}

SpectralReal DynamicAdversary::objective_value(
    const EvolutionResult& evolution, const std::string& objective) {
    if (objective == "critical-integral") {
        return evolution.integral_critical;
    }
    if (objective == "critical-local-integral") {
        return evolution.integral_local_critical;
    }
    if (objective == "critical-local-increase") {
        return evolution.final_local_critical_integrand -
               evolution.initial_local_critical_integrand;
    }
    if (objective == "critical-nonlocal-integral") {
        return evolution.integral_nonlocal_critical;
    }
    if (objective == "critical-near-nonlocal-integral") {
        return evolution.integral_near_nonlocal_critical;
    }
    if (objective == "critical-far-nonlocal-integral") {
        return evolution.integral_far_nonlocal_critical;
    }
    if (objective == "critical-gap-tail-integral") {
        return evolution.integral_selected_gap_tail_critical;
    }
    if (objective == "max-q") {
        return evolution.maximum_energy_level_quantity;
    }
    if (objective == "terminal-q") {
        return evolution.final_energy_level_quantity;
    }
    if (objective == "q-gain") {
        if (!(evolution.initial_energy_level_quantity > 1e-30L) ||
            !(evolution.final_energy_level_quantity > 1e-30L)) {
            return -std::numeric_limits<SpectralReal>::infinity();
        }
        return std::log(evolution.final_energy_level_quantity /
                        evolution.initial_energy_level_quantity);
    }
    if (objective == "q-increase") {
        return evolution.final_energy_level_quantity -
               evolution.initial_energy_level_quantity;
    }
    throw std::invalid_argument("unknown dynamic objective: " + objective);
}

DynamicAdversaryResult DynamicAdversary::optimize(
    const SpectralState& primary_start,
    const SpectralState* secondary_start,
    const DynamicAdversaryOptions& options,
    bool refine_result) const {
    if (options.generations < 0) {
        throw std::invalid_argument("--dynamic-generations cannot be negative");
    }
    std::mt19937_64 generator(
        options.seed ^ UINT64_C(0xd1b54a32d192ed03) ^
        static_cast<std::uint64_t>(SpectralStateOps::cutoff(primary_start)));
    DynamicAdversaryResult result;
    const InitialSobolevConstraint sobolev(
        options.sobolev_order, options.sobolev_cap);
    const bool collect_search_partition =
        collects_partition(options.objective);
    result.state = primary_start;
    result.initial_objective = trajectory_.evaluate_static(result.state);
    result.evolution = trajectory_.evolve(
        result.state, options.viscosity, options.final_time,
        options.time_step, collect_search_partition,
        options.minimum_dyadic_gap);
    ++result.evaluations;
    bool result_admissible = sobolev.admissible(result.state);

    if (secondary_start != nullptr) {
        SpectralState secondary = move_to_cutoff(
            *secondary_start, SpectralStateOps::cutoff(primary_start),
            generator);
        const EvolutionResult secondary_evolution = trajectory_.evolve(
            secondary, options.viscosity, options.final_time,
            options.time_step, collect_search_partition,
            options.minimum_dyadic_gap);
        ++result.evaluations;
        const bool secondary_admissible = sobolev.admissible(secondary);
        if (secondary_admissible &&
            (!result_admissible ||
             objective_value(secondary_evolution, options.objective) >
                 objective_value(result.evolution, options.objective))) {
            result.state = std::move(secondary);
            result.initial_objective =
                trajectory_.evaluate_static(result.state);
            result.evolution = secondary_evolution;
            result_admissible = true;
        }
    }
    if (!result_admissible) {
        throw std::invalid_argument(
            "no dynamic start satisfies the configured Sobolev cap");
    }
    result.search_initial_objective =
        objective_value(result.evolution, options.objective);

    const bool use_mutations = options.optimizer == "mutate" ||
                               options.optimizer == "hybrid";
    const bool use_gradient = options.optimizer == "gradient" ||
                              options.optimizer == "hybrid";
    for (int generation = 0;
         use_mutations && generation < options.generations;
         ++generation) {
        const SpectralReal progress =
            options.generations > 1
                ? static_cast<SpectralReal>(generation) /
                      static_cast<SpectralReal>(options.generations - 1)
                : 0.0L;
        const SpectralReal radius =
            options.mutation * (0.5L - 0.4L * progress);
        SpectralState candidate = SpectralStateFactory::mutate(
            result.state, radius, generator, generation % 4 != 0);
        if (!sobolev.admissible(candidate)) {
            continue;
        }
        const EvolutionResult candidate_evolution = trajectory_.evolve(
            candidate, options.viscosity, options.final_time,
            options.time_step, collect_search_partition,
            options.minimum_dyadic_gap);
        ++result.evaluations;
        if (candidate_evolution.finite &&
            objective_value(candidate_evolution, options.objective) >
                objective_value(result.evolution, options.objective)) {
            result.state = std::move(candidate);
            result.initial_objective =
                trajectory_.evaluate_static(result.state);
            result.evolution = candidate_evolution;
            ++result.accepted_mutations;
        }
    }
    if (use_gradient && options.generations > 0) {
        const int trajectory_steps = std::max(
            1, static_cast<int>(
                   std::ceil(options.final_time / options.time_step)));
        GradientSearchOptions gradient_options;
        gradient_options.iterations = options.generations;
        gradient_options.line_search_steps = 16;
        gradient_options.trajectory_steps = trajectory_steps;
        gradient_options.viscosity = options.viscosity;
        gradient_options.time_step =
            options.final_time / static_cast<SpectralReal>(trajectory_steps);
        gradient_options.initial_step = options.mutation;
        gradient_options.objective = options.objective;
        gradient_options.method = options.gradient_method;
        gradient_options.sobolev_order = options.sobolev_order;
        gradient_options.sobolev_cap = options.sobolev_cap;
        gradient_options.minimum_dyadic_gap = options.minimum_dyadic_gap;
        const GradientSearchResult gradient =
            gradient_.maximize_q(result.state, gradient_options);
        result.state = gradient.state;
        result.initial_objective = trajectory_.evaluate_static(result.state);
        result.evolution = trajectory_.evolve(
            result.state, options.viscosity, options.final_time,
            options.time_step, collect_search_partition,
            options.minimum_dyadic_gap);
        result.evaluations += gradient.trajectory_evaluations + 1;
        result.accepted_gradient_steps = gradient.accepted_steps;
        result.gradient_trace = gradient.trace;
    }
    result.search_final_objective =
        objective_value(result.evolution, options.objective);
    if (refine_result) {
        refine(result, options);
    }
    return result;
}

void DynamicAdversary::refine(
    DynamicAdversaryResult& result,
    const DynamicAdversaryOptions& options) const {
    result.refined_evolution = trajectory_.evolve(
        result.state, options.viscosity, options.final_time,
        0.5L * options.time_step, true, options.minimum_dyadic_gap);
    const SpectralReal refined_objective = objective_value(
        result.refined_evolution, options.objective);
    const SpectralReal coarse_objective = objective_value(
        result.evolution, options.objective);
    const SpectralReal difference =
        std::abs(refined_objective - coarse_objective);
    result.time_step_relative_error =
        refined_objective != 0.0L
            ? difference / std::abs(refined_objective)
            : (difference == 0.0L
                   ? 0.0L
                   : std::numeric_limits<SpectralReal>::infinity());
}

DynamicAdversaryEnsemble::DynamicAdversaryEnsemble(
    std::string backend, int workers)
    : backend_(std::move(backend)), workers_(workers) {
    if (workers_ < 1 || workers_ > 256) {
        throw std::invalid_argument(
            "dynamic adversary workers must be between 1 and 256");
    }
}

DynamicAdversaryResult DynamicAdversaryEnsemble::optimize(
    const SpectralState& primary_start,
    const SpectralState* secondary_start,
    const DynamicAdversaryOptions& options,
    int restarts) const {
    if (restarts < 1 || restarts > 10000) {
        throw std::invalid_argument(
            "dynamic adversary restarts must be between 1 and 10000");
    }
    if (restarts == 1) {
        DynamicAdversary adversary(backend_, workers_);
        DynamicAdversaryResult result = adversary.optimize(
            primary_start, secondary_start, options, true);
        result.restart_objectives = {result.search_final_objective};
        return result;
    }

    const int cutoff = SpectralStateOps::cutoff(primary_start);
    const SpectralReal target_energy =
        SpectralStateOps::energy(primary_start);
    const InitialSobolevConstraint sobolev(
        options.sobolev_order, options.sobolev_cap);
    std::mt19937_64 layout_generator(options.seed ^ UINT64_C(0x6a09e667f3bcc909));
    SpectralState warm = secondary_start != nullptr
        ? move_to_cutoff(*secondary_start, cutoff, layout_generator)
        : primary_start;
    sobolev.retract(warm, target_energy);

    std::vector<SpectralState> starts(static_cast<std::size_t>(restarts));
    starts.front() = primary_start;
    for (int restart = 1; restart < restarts; ++restart) {
        std::mt19937_64 generator(
            options.seed + static_cast<std::uint64_t>(restart) *
                               UINT64_C(0x94d049bb133111eb));
        SpectralState start;
        if (restart % 2 == 1) {
            const SpectralReal radius = options.mutation *
                (2.0L + static_cast<SpectralReal>(restart) /
                            static_cast<SpectralReal>(restarts));
            start = SpectralStateFactory::mutate(
                warm, radius, generator, restart % 4 != 1);
        } else {
            start = SpectralStateFactory::random(cutoff, generator);
        }
        sobolev.retract(start, target_energy);
        starts[static_cast<std::size_t>(restart)] = std::move(start);
    }

    // Populate the process-wide interaction cache before worker threads only
    // read it. Spectral states with the same cutoff share this layout.
    static_cast<void>(SpectralStateOps::interactions(primary_start));
    std::vector<DynamicAdversaryResult> results(
        static_cast<std::size_t>(restarts));
    const ParallelExecutor executor(workers_);
    executor.for_each(results.size(), [&](std::size_t restart) {
        DynamicAdversaryOptions local_options = options;
        local_options.seed += static_cast<std::uint64_t>(restart) *
                              UINT64_C(0xbf58476d1ce4e5b9);
        DynamicAdversary adversary(backend_, 1);
        results[restart] = adversary.optimize(
            starts[restart], restart == 0 ? secondary_start : nullptr,
            local_options, false);
    });

    std::size_t winner = 0;
    int total_evaluations = 0;
    std::vector<SpectralReal> restart_objectives(results.size());
    for (std::size_t restart = 0; restart < results.size(); ++restart) {
        total_evaluations += results[restart].evaluations;
        restart_objectives[restart] =
            results[restart].search_final_objective;
        if (results[restart].search_final_objective >
            results[winner].search_final_objective) {
            winner = restart;
        }
    }
    DynamicAdversaryResult best = std::move(results[winner]);
    best.evaluations = total_evaluations;
    best.winning_restart = static_cast<int>(winner);
    best.restart_objectives = std::move(restart_objectives);
    DynamicAdversary refiner(backend_, workers_);
    refiner.refine(best, options);
    return best;
}

}  // namespace lemma
