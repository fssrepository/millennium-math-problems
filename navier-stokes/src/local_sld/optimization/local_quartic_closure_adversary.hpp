#pragma once

#include "gradient_adversary.hpp"
#include "local_quartic_closure_objective.hpp"
#include "local_sld_block_objective.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace lemma {

struct LocalQuarticClosureAdversaryOptions {
    int minimum_cutoff = 2;
    int maximum_cutoff = 6;
    int restarts = 12;
    int workers = 12;
    int iterations = 8;
    int line_search_steps = 14;
    int lbfgs_history = 8;
    int trajectory_steps = 8;
    SpectralReal initial_step = 0.1L;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.002L;
    SpectralReal absorption_theta = 0.9L;
    int shape_power = 3;
    int sobolev_order = 0;
    SpectralReal sobolev_cap = 0.0L;
    std::uint64_t seed = 20260801;
    std::string method = "lbfgs";
    std::string backend = "direct";
    std::string objective = "sld-ratio";
    std::string selection = "local";
    std::string initial_profile = "mixed";
    std::string certificate_path;
    std::string state_directory;
    std::string warm_state_path;
};

struct LocalQuarticClosureRestartResult {
    SpectralState state;
    LocalQuarticClosureObjectiveValue value;
    LocalSldBlockObjectiveValue common_block_value;
    bool common_block_objective = false;
    SpectralReal refined_objective = 0.0L;
    SpectralReal time_step_relative_error = 0.0L;
    SpectralReal frozen_initial_frequency = 0.0L;
    SpectralReal frozen_initial_ep_shift = 0.0L;
    int objective_step = 0;
    int refined_objective_step = 0;
    SpectralReal initial_objective = 0.0L;
    SpectralReal objective = 0.0L;
    SpectralReal initial_constant_ratio = 0.0L;
    SpectralReal remainder_envelope_ratio = 0.0L;
    SpectralReal remainder_absorption_ratio = 0.0L;
    SpectralReal shape_power_absolute_product = 0.0L;
    SpectralReal shape_power_normalized_stretching = 0.0L;
    SpectralReal projective_coherence_ratio = 0.0L;
    SpectralReal projective_coherence_amplification = 0.0L;
    std::size_t projective_coherence_shape_count = 0;
    SpectralReal final_projected_gradient_norm = 0.0L;
    SpectralReal sobolev_value = 0.0L;
    std::uint64_t seed = 0;
    int restart = 0;
    int accepted_steps = 0;
    int evaluations = 0;
    bool warm_continuation = false;
};

struct LocalQuarticClosureCutoffResult {
    int cutoff = 0;
    LocalQuarticClosureRestartResult winner;
    SpectralReal warm_lift_constant_ratio = 0.0L;
    SpectralReal warm_lift_objective = 0.0L;
    SpectralReal improvement_factor = 0.0L;
    SpectralReal objective_gain = 0.0L;
    SpectralReal projection_residual = 0.0L;
    std::string state_path;
    std::vector<SpectralReal> restart_constant_ratios;
    std::vector<SpectralReal> restart_objectives;
};

struct LocalQuarticClosureAdversaryReport {
    int workers = 0;
    int restarts = 0;
    int iterations = 0;
    int trajectory_steps = 0;
    std::string objective;
    std::string backend;
    std::string initial_profile;
    SpectralReal viscosity = 0.0L;
    SpectralReal time_step = 0.0L;
    SpectralReal absorption_theta = 0.0L;
    int shape_power = 0;
    int sobolev_order = 0;
    SpectralReal sobolev_cap = 0.0L;
    SpectralReal fitted_cutoff_slope = 0.0L;
    SpectralReal maximum_constant_ratio = 0.0L;
    SpectralReal maximum_objective = 0.0L;
    bool finite_search_is_not_a_proof = true;
    bool candidate_lemma_proved = false;
    std::vector<LocalQuarticClosureCutoffResult> rows;
};

class LocalQuarticClosureAdversary {
public:
    [[nodiscard]] static LocalQuarticClosureRestartResult maximize(
        const SpectralState& initial,
        const LocalQuarticClosureAdversaryOptions& options,
        int restart, std::uint64_t seed, bool warm_continuation);
};

class LocalQuarticClosureEnsemble {
public:
    [[nodiscard]] static LocalQuarticClosureAdversaryReport scan(
        const LocalQuarticClosureAdversaryOptions& options);
};

}  // namespace lemma
