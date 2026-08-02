#include "local_quartic_closure_adversary.hpp"

#include "initial_sobolev_constraint.hpp"
#include "local_signature_state_factory.hpp"
#include "local_sld_trajectory_adjoint.hpp"
#include "local_sld_remainder_envelope_objective.hpp"
#include "local_sld_remainder_absorption_objective.hpp"
#include "local_sld_shape_power_objective.hpp"
#include "local_sld_projective_coherence_objective.hpp"
#include "local_sld_projective_stretching_objective.hpp"
#include "local_sld_projective_cross_power_objective.hpp"
#include "local_sld_projective_open_power_objective.hpp"
#include "local_sld_projective_height_stretching_objective.hpp"
#include "local_sld_projective_height_power_objective.hpp"
#include "local_sld_projective_height_outer_power_objective.hpp"
#include "local_sld_projective_height_envelope_objective.hpp"
#include "local_sld_projective_height_commutator_ratio_objective.hpp"
#include "local_sld_projective_height_dynamic_ratio_objective.hpp"
#include "local_sld_projective_normalization_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "parallel_executor.hpp"
#include "spectral_adjoint.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_objective.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace lemma {
namespace {

TriadSelection closure_selection(const std::string& name) {
    return LocalSldTriadSelection::parse(name);
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
    const InitialSobolevConstraint& sobolev,
    const std::string& initial_profile,
    bool preserve_warm_layout) {
    std::mt19937_64 generator(seed);
    SpectralState state;
    if (previous_winner != nullptr && restart == 0) {
        state = preserve_warm_layout &&
                SpectralStateOps::cutoff(*previous_winner) == cutoff
            ? *previous_winner
            : SpectralStateFactory::lift(
                  *previous_winner, cutoff, generator);
    } else if (initial_profile != "mixed") {
        state = LocalSignatureStateFactory::make(
            cutoff,
            LocalSignatureStateFactory::parse(initial_profile),
            seed);
    } else if (previous_winner != nullptr && restart % 3 == 0) {
        state = preserve_warm_layout &&
                SpectralStateOps::cutoff(*previous_winner) == cutoff
            ? *previous_winner
            : SpectralStateFactory::lift(
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

void write_restart_checkpoint(
    const LocalQuarticClosureAdversaryOptions& options,
    int cutoff,
    const LocalQuarticClosureRestartResult& result) {
    if (options.state_directory.empty()) {
        return;
    }
    const std::filesystem::path directory =
        std::filesystem::path(options.state_directory) /
        "restarts" / ("K" + std::to_string(cutoff));
    std::filesystem::create_directories(directory);
    std::ostringstream filename;
    filename << 'R' << std::setw(3) << std::setfill('0')
             << result.restart << ".tsv";
    std::ostringstream metadata;
    metadata << std::setprecision(18)
             << "completed local quartic closure restart checkpoint; "
             << "cutoff=" << cutoff
             << "; restart=" << result.restart
             << "; seed=" << result.seed
             << "; objective=" << options.objective
             << "; selection=" << options.selection
             << "; optimized_objective="
             << static_cast<double>(result.objective)
             << "; refined_objective="
             << static_cast<double>(result.refined_objective)
             << "; accepted_steps=" << result.accepted_steps
             << "; candidate_lemma_proved=false";
    SpectralStateWriter::write_tsv(
        (directory / filename.str()).string(),
        result.state, metadata.str());
}

}  // namespace

LocalQuarticClosureRestartResult
LocalQuarticClosureAdversary::maximize(
    const SpectralState& initial,
    const LocalQuarticClosureAdversaryOptions& options,
    int restart, std::uint64_t seed, bool warm_continuation) {
    SpectralGalerkin galerkin;
    galerkin.configure(
        options.backend,
        options.restarts == 1 ? options.workers : 1);
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
    search.absorption_theta = options.absorption_theta;
    search.shape_power = options.shape_power;
    search.projective_core_maximum_height =
        options.projective_core_maximum_height;
    search.objective = options.objective == "closure-ratio"
        ? "local-closure-ratio"
        : (options.objective == "lqc3-ratio"
               ? "local-lqc3-ratio"
        : (options.objective == "signed-lqc3-ratio"
               ? "local-signed-lqc3-ratio"
        : (options.objective == "remainder-envelope-ratio"
               ? "local-remainder-envelope-ratio"
        : (options.objective == "remainder-absorption-ratio"
               ? "local-remainder-absorption-ratio"
        : (options.objective == "shape-power-ratio"
               ? "local-shape-power-ratio"
        : (options.objective == "projective-coherence-ratio"
               ? "local-projective-coherence-ratio"
        : (options.objective == "projective-stretching-ratio"
               ? "local-projective-stretching-ratio"
        : (options.objective == "projective-cross-power-ratio"
               ? "local-projective-cross-power-ratio"
        : (options.objective == "projective-open-power-ratio"
               ? "local-projective-open-power-ratio"
        : (options.objective == "projective-height-stretching-ratio"
               ? "local-projective-height-stretching-ratio"
        : (options.objective == "projective-height-power-ratio"
               ? "local-projective-height-power-ratio"
        : (options.objective == "projective-height-outer-power-ratio"
               ? "local-projective-height-outer-power-ratio"
        : (options.objective == "projective-height-envelope-ratio"
               ? "local-projective-height-envelope-ratio"
        : (options.objective ==
               "projective-height-commutator-envelope-ratio"
               ? "local-projective-height-commutator-envelope-ratio"
        : (options.objective ==
               "projective-height-commutator-coercivity-ratio"
               ? "local-projective-height-commutator-coercivity-ratio"
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
                                           : "local-sld-ratio"))))))))))))))))))));
    search.method = options.method;
    if (options.objective ==
        "projective-height-dynamic-coercivity-ratio") {
        search.objective =
            "local-projective-height-dynamic-coercivity-ratio";
    }
    if (options.objective ==
        "projective-height-dynamic-envelope-ratio") {
        search.objective =
            "local-projective-height-dynamic-envelope-ratio";
    }
    if (options.objective ==
        "projective-palinstrophy-normalization-ratio") {
        search.objective =
            "local-projective-palinstrophy-normalization-ratio";
    }
    search.lbfgs_history = options.lbfgs_history;
    search.sobolev_order = options.sobolev_order;
    search.sobolev_cap = options.sobolev_cap;
    search.closure_selection = closure_selection(options.selection);
    search.objective_threads = options.restarts == 1
        ? options.workers : 1;
    const GradientSearchResult optimized = adversary.maximize_q(
        initial, search);
    LocalQuarticClosureRestartResult result;
    result.state = optimized.state;
    if (options.lean_diagnostics) {
        result.value.finite = std::isfinite(optimized.objective);
        if (options.objective ==
            "projective-height-envelope-ratio") {
            result.projective_height_component_envelope_absolute =
                std::sqrt(std::max(0.0L, optimized.objective));
        }
        if (options.objective ==
            "projective-height-commutator-envelope-ratio") {
            result.projective_height_commutator_envelope_absolute =
                std::sqrt(std::max(0.0L, optimized.objective));
        }
        if (options.objective ==
            "projective-height-dynamic-envelope-ratio") {
            result.projective_height_dynamic_envelope_absolute =
                std::sqrt(std::max(0.0L, optimized.objective));
        }
        if (options.objective ==
            "projective-height-commutator-coercivity-ratio") {
            result.projective_height_commutator_coercivity_ratio =
                std::sqrt(std::max(0.0L, optimized.objective));
        }
        if (options.objective ==
            "projective-height-dynamic-coercivity-ratio") {
            result.projective_height_dynamic_coercivity_ratio =
                std::sqrt(std::max(0.0L, optimized.objective));
        }
        if (options.objective ==
            "projective-palinstrophy-normalization-ratio") {
            result.projective_palinstrophy_normalization_power_one =
                std::sqrt(std::max(0.0L, optimized.objective));
        }
        result.initial_objective = optimized.initial_objective;
        result.objective = optimized.objective;
        result.objective_step = optimized.objective_step;
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
    const LocalQuarticClosureObjective closure(
        dynamics, search.closure_selection);
    const LocalQuarticClosureObjectiveValue initial_value =
        closure.evaluate(initial);
    result.value = closure.evaluate(result.state);
    result.remainder_envelope_ratio =
        LocalSldRemainderEnvelopeObjective(
            dynamics, search.closure_selection)
            .evaluate(result.state).target_ratio;
    result.remainder_absorption_ratio =
        LocalSldRemainderAbsorptionObjective(
            dynamics, options.absorption_theta,
            search.closure_selection)
            .evaluate(result.state).absorption_ratio;
    const LocalSldShapePowerObjectiveValue shape_power_value =
        LocalSldShapePowerObjective(
            dynamics, search.closure_selection,
            options.shape_power).evaluate(result.state);
    result.shape_power_absolute_product =
        shape_power_value.absolute_power_product;
    result.shape_power_normalized_stretching =
        shape_power_value.normalized_stretching;
    const LocalSldProjectiveCoherenceObjectiveValue coherence_value =
        LocalSldProjectiveCoherenceObjective(
            dynamics, search.closure_selection).evaluate(result.state);
    result.projective_coherence_ratio = coherence_value.synthesis_ratio;
    result.projective_coherence_amplification =
        coherence_value.synthesis_amplification;
    result.projective_coherence_shape_count =
        coherence_value.projective_shape_count;
    const LocalSldProjectiveStretchingObjectiveValue stretching_value =
        LocalSldProjectiveStretchingObjective(
            dynamics, search.closure_selection).evaluate(result.state);
    result.projective_stretching_ratio =
        stretching_value.stretching_aware_synthesis_ratio;
    result.projective_stretching_alignment_squared =
        stretching_value.stretching_alignment_squared;
    result.projective_stretching_reconstruction_error =
        stretching_value.product_reconstruction_error;
    const LocalSldProjectiveCrossPowerObjectiveValue cross_power_value =
        LocalSldProjectiveCrossPowerObjective(
            dynamics, search.closure_selection,
            search.objective_threads).evaluate(result.state);
    result.projective_cross_power_absolute =
        cross_power_value.absolute_cross_power_one;
    result.projective_cross_bracket = cross_power_value.cross_bracket;
    result.projective_diagonal_bracket =
        cross_power_value.diagonal_bracket;
    const LocalSldProjectiveOpenPowerObjectiveValue open_power_value =
        LocalSldProjectiveOpenPowerObjective(
            dynamics, search.closure_selection,
            options.projective_core_maximum_height,
            search.objective_threads).evaluate(result.state);
    result.projective_open_power_absolute =
        open_power_value.absolute_open_power_one;
    result.projective_open_bracket = open_power_value.open_bracket;
    result.projective_fixed_core_bracket =
        open_power_value.fixed_core_bracket;
    const LocalSldProjectiveHeightStretchingObjectiveValue
        height_stretching_value =
            LocalSldProjectiveHeightStretchingObjective(
                dynamics, search.closure_selection,
                options.projective_core_maximum_height,
                search.objective_threads).evaluate(result.state);
    result.projective_height_stretching_ratio =
        height_stretching_value.stretching_aware_h1_ratio;
    result.projective_height_h1_synthesis_ratio =
        height_stretching_value.h1_synthesis_ratio;
    result.projective_height_stretching_alignment_squared =
        height_stretching_value.stretching_h1_alignment_squared;
    result.projective_height_shape_count =
        height_stretching_value.shell_shape_count;
    const LocalSldProjectiveHeightPowerObjectiveValue
        height_power_value = LocalSldProjectiveHeightPowerObjective(
            dynamics, search.closure_selection,
            options.projective_core_maximum_height,
            search.objective_threads).evaluate(result.state);
    result.projective_height_power_absolute =
        height_power_value.absolute_shell_power_one;
    result.projective_height_internal_bracket =
        height_power_value.shell_internal_bracket;
    const LocalSldProjectiveHeightOuterPowerObjectiveValue
        height_outer_power_value =
            LocalSldProjectiveHeightOuterPowerObjective(
                dynamics, search.closure_selection,
                search.objective_threads).evaluate(result.state);
    result.projective_height_outer_power_absolute =
        height_outer_power_value.absolute_outer_power_one;
    result.projective_height_outer_h1_sum =
        height_outer_power_value.diagonal_outer_h1_sum;
    result.projective_height_active_shell_count =
        height_outer_power_value.active_height_shell_count;
    const LocalSldProjectiveHeightEnvelopeObjectiveValue
        height_envelope_value =
            LocalSldProjectiveHeightEnvelopeObjective(
                dynamics, search.closure_selection,
                search.objective_threads).evaluate(result.state);
    result.projective_height_component_envelope_absolute =
        height_envelope_value.absolute_component_power_one_envelope;
    result.projective_height_commutator_envelope_absolute =
        LocalSldProjectiveHeightEnvelopeObjective(
            dynamics, search.closure_selection,
            search.objective_threads, true)
            .evaluate(result.state)
            .absolute_component_power_one_envelope;
    result.projective_height_dynamic_envelope_absolute =
        LocalSldProjectiveHeightEnvelopeObjective(
            dynamics, search.closure_selection,
            search.objective_threads, true, true)
            .evaluate(result.state)
            .absolute_component_power_one_envelope;
    result.projective_height_commutator_coercivity_ratio =
        LocalSldProjectiveHeightCommutatorRatioObjective(
            dynamics, search.closure_selection,
            search.objective_threads)
            .evaluate(result.state).coercivity_ratio;
    result.projective_height_dynamic_coercivity_ratio =
        LocalSldProjectiveHeightDynamicRatioObjective(
            dynamics, search.closure_selection,
            search.objective_threads)
            .evaluate(result.state).coercivity_ratio;
    result.projective_palinstrophy_normalization_power_one =
        LocalSldProjectiveNormalizationObjective(
            dynamics, search.closure_selection)
            .evaluate(result.state)
            .palinstrophy_normalization_power_one;
    result.projective_height_component_bracket_envelope =
        height_envelope_value.absolute_component_bracket_envelope;
    result.projective_height_pair_count =
        height_envelope_value.height_pair_count;
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
        options.maximum_cutoff > 16 || options.restarts < 1 ||
        options.restarts > 1000 || options.workers < 1 ||
        options.workers > 256 || options.iterations < 0 ||
        options.line_search_steps < 1 ||
        !(options.initial_step > 0.0L) ||
        !(options.absorption_theta >= 0.0L) ||
        !(options.absorption_theta <= 1.0L) ||
        !std::isfinite(options.absorption_theta) ||
        options.shape_power < 0 || options.shape_power > 3 ||
        options.projective_core_maximum_height < 1 ||
        options.projective_core_maximum_height > 256 ||
        (options.method != "steepest" && options.method != "lbfgs") ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft") ||
        (options.objective != "closure-ratio" &&
         options.objective != "lqc3-ratio" &&
         options.objective != "signed-lqc3-ratio" &&
         options.objective != "remainder-envelope-ratio" &&
         options.objective != "remainder-absorption-ratio" &&
         options.objective != "shape-power-ratio" &&
         options.objective != "projective-coherence-ratio" &&
         options.objective != "projective-stretching-ratio" &&
         options.objective != "projective-cross-power-ratio" &&
         options.objective != "projective-open-power-ratio" &&
         options.objective != "projective-height-stretching-ratio" &&
         options.objective != "projective-height-power-ratio" &&
         options.objective != "projective-height-outer-power-ratio" &&
         options.objective != "projective-height-envelope-ratio" &&
         options.objective !=
             "projective-height-commutator-envelope-ratio" &&
         options.objective !=
             "projective-height-dynamic-envelope-ratio" &&
         options.objective !=
             "projective-height-commutator-coercivity-ratio" &&
         options.objective !=
             "projective-height-dynamic-coercivity-ratio" &&
         options.objective !=
             "projective-palinstrophy-normalization-ratio" &&
         options.objective != "signed-closure-ratio" &&
         options.objective != "sld-ratio" &&
         options.objective != "block-ratio" &&
         options.objective != "mixed-ratio" &&
         options.objective != "terminal-sld-ratio" &&
         options.objective != "maximum-sld-ratio") ||
        !LocalSldTriadSelection::supports(options.selection) ||
        (options.initial_profile != "mixed" &&
         options.initial_profile != "decaying" &&
         options.initial_profile != "flat" &&
         options.initial_profile != "outer-half-flat") ||
        (is_common_block_objective(options.objective) &&
         (options.selection == "local" ||
          options.selection == "remainder-without-123" ||
          options.selection ==
              "double-triple-remainder-without-123")) ||
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
    report.initial_profile = options.initial_profile;
    report.viscosity = options.viscosity;
    report.time_step = options.time_step;
    report.absorption_theta = options.absorption_theta;
    report.shape_power = options.shape_power;
    report.projective_core_maximum_height =
        options.projective_core_maximum_height;
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
                sobolev, options.initial_profile,
                options.preserve_warm_layout);
        }

        LocalQuarticClosureCutoffResult row;
        row.cutoff = cutoff;
        if (has_previous_winner) {
            SpectralGalerkin galerkin;
            galerkin.configure(options.backend, 1);
            const SpectralDynamics dynamics(galerkin);
            if (options.lean_diagnostics &&
                (options.objective ==
                     "projective-height-envelope-ratio" ||
                 options.objective ==
                     "projective-height-commutator-envelope-ratio" ||
                 options.objective ==
                     "projective-height-dynamic-envelope-ratio" ||
                 options.objective ==
                     "projective-height-commutator-coercivity-ratio" ||
                 options.objective ==
                     "projective-height-dynamic-coercivity-ratio" ||
                 options.objective ==
                     "projective-palinstrophy-normalization-ratio")) {
                if (options.objective ==
                    "projective-palinstrophy-normalization-ratio") {
                    row.warm_lift_objective =
                        LocalSldProjectiveNormalizationObjective(
                            dynamics,
                            closure_selection(options.selection))
                            .evaluate(starts.front())
                            .squared_palinstrophy_normalization_power_one;
                } else if (options.objective ==
                    "projective-height-commutator-coercivity-ratio") {
                    row.warm_lift_objective =
                        LocalSldProjectiveHeightCommutatorRatioObjective(
                            dynamics,
                            closure_selection(options.selection),
                            options.workers)
                            .evaluate(starts.front())
                            .squared_coercivity_ratio;
                } else if (options.objective ==
                           "projective-height-dynamic-coercivity-ratio") {
                    row.warm_lift_objective =
                        LocalSldProjectiveHeightDynamicRatioObjective(
                            dynamics,
                            closure_selection(options.selection),
                            options.workers)
                            .evaluate(starts.front())
                            .squared_coercivity_ratio;
                } else {
                    row.warm_lift_objective =
                        LocalSldProjectiveHeightEnvelopeObjective(
                            dynamics,
                            closure_selection(options.selection),
                            options.workers,
                            options.objective ==
                                "projective-height-commutator-envelope-ratio" ||
                            options.objective ==
                                "projective-height-dynamic-envelope-ratio",
                            options.objective ==
                                "projective-height-dynamic-envelope-ratio")
                            .evaluate(starts.front())
                            .squared_component_power_one_envelope;
                }
            } else {
            const LocalQuarticClosureObjectiveValue warm_value =
                LocalQuarticClosureObjective(
                    dynamics, closure_selection(options.selection))
                    .evaluate(starts.front());
            row.warm_lift_constant_ratio = warm_value.constant_ratio;
            row.warm_lift_objective = options.objective == "closure-ratio"
                ? warm_value.squared_constant_ratio
                : (options.objective == "lqc3-ratio"
                       ? warm_value.squared_lqc3_target_ratio
                : (options.objective == "signed-lqc3-ratio"
                       ? warm_value.signed_lqc3_target_ratio
                : (options.objective == "remainder-envelope-ratio"
                       ? LocalSldRemainderEnvelopeObjective(
                             dynamics,
                             closure_selection(options.selection))
                             .evaluate(starts.front()).target_ratio
                : (options.objective == "remainder-absorption-ratio"
                       ? LocalSldRemainderAbsorptionObjective(
                             dynamics, options.absorption_theta,
                             closure_selection(options.selection))
                             .evaluate(starts.front()).absorption_ratio
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
                              : warm_value.signed_local_sld_ratio))))))));
            if (options.objective == "projective-open-power-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveOpenPowerObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.projective_core_maximum_height)
                        .evaluate(starts.front())
                        .squared_open_power_one;
            }
            if (options.objective ==
                "projective-height-stretching-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightStretchingObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.projective_core_maximum_height)
                        .evaluate(starts.front())
                        .stretching_aware_h1_ratio;
            }
            if (options.objective ==
                "projective-height-power-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightPowerObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.projective_core_maximum_height)
                        .evaluate(starts.front())
                        .squared_shell_power_one;
            }
            if (options.objective ==
                "projective-height-outer-power-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightOuterPowerObjective(
                        dynamics,
                        closure_selection(options.selection))
                        .evaluate(starts.front())
                        .squared_outer_power_one;
            }
            if (options.objective ==
                "projective-height-envelope-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightEnvelopeObjective(
                        dynamics,
                        closure_selection(options.selection))
                        .evaluate(starts.front())
                        .squared_component_power_one_envelope;
            }
            if (options.objective ==
                "projective-height-commutator-envelope-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightEnvelopeObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.workers, true)
                        .evaluate(starts.front())
                        .squared_component_power_one_envelope;
            }
            if (options.objective ==
                "projective-height-dynamic-envelope-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightEnvelopeObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.workers, true, true)
                        .evaluate(starts.front())
                        .squared_component_power_one_envelope;
            }
            if (options.objective ==
                "projective-height-commutator-coercivity-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightCommutatorRatioObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.workers)
                        .evaluate(starts.front())
                        .squared_coercivity_ratio;
            }
            if (options.objective ==
                "projective-height-dynamic-coercivity-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveHeightDynamicRatioObjective(
                        dynamics,
                        closure_selection(options.selection),
                        options.workers)
                        .evaluate(starts.front())
                        .squared_coercivity_ratio;
            }
            if (options.objective ==
                "projective-palinstrophy-normalization-ratio") {
                row.warm_lift_objective =
                    LocalSldProjectiveNormalizationObjective(
                        dynamics,
                        closure_selection(options.selection))
                        .evaluate(starts.front())
                        .squared_palinstrophy_normalization_power_one;
            }
            }
        }
        std::vector<LocalQuarticClosureRestartResult> results(
            static_cast<std::size_t>(options.restarts));
        executor.for_each(results.size(), [&](std::size_t restart) {
            results[restart] = LocalQuarticClosureAdversary::maximize(
                starts[restart], options, static_cast<int>(restart),
                seeds[restart], has_previous_winner && restart == 0);
            write_restart_checkpoint(
                options, cutoff, results[restart]);
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
