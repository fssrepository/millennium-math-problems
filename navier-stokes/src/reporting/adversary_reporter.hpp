#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct AdversaryGradientTracePoint {
    int iteration = 0;
    long double objective_before = 0.0L;
    long double objective_after = 0.0L;
    long double projected_gradient_norm = 0.0L;
    long double accepted_step = 0.0L;
    long double sobolev_value = 0.0L;
    int line_search_evaluations = 0;
    bool used_steepest_fallback = false;
    bool accepted = false;
};

struct AdversaryReportRow {
    int cutoff = 0;
    int modes = 0;
    int evaluations = 0;
    int accepted_mutations = 0;
    long double energy = 0.0L;
    long double enstrophy = 0.0L;
    long double palinstrophy = 0.0L;
    long double vortex_stretching = 0.0L;
    long double depletion = 0.0L;
    long double q = 0.0L;
    long double critical_integrand = 0.0L;
    int dynamic_steps = 0;
    long double dynamic_integral = 0.0L;
    long double dynamic_coarse_integral = 0.0L;
    long double dynamic_local_integral = 0.0L;
    long double dynamic_nonlocal_integral = 0.0L;
    long double dynamic_near_nonlocal_integral = 0.0L;
    long double dynamic_far_nonlocal_integral = 0.0L;
    long double dynamic_selected_gap_tail_integral = 0.0L;
    long double dynamic_dt_relative_error = 0.0L;
    long double dynamic_search_initial_objective = 0.0L;
    long double dynamic_search_final_objective = 0.0L;
    long double dynamic_initial_local_critical_density = 0.0L;
    long double dynamic_final_local_critical_density = 0.0L;
    long double dynamic_initial_enstrophy = 0.0L;
    long double dynamic_initial_palinstrophy = 0.0L;
    long double dynamic_initial_ep_shift = 0.0L;
    long double dynamic_local_critical_log_gain = 0.0L;
    long double dynamic_local_log_gain_rate_ratio = 0.0L;
    long double dynamic_ep_shifted_local_log_gain = 0.0L;
    long double dynamic_ep_shifted_log_gain_rate_ratio = 0.0L;
    long double dynamic_maximum_q = 0.0L;
    long double dynamic_initial_q = 0.0L;
    long double dynamic_final_q = 0.0L;
    long double dynamic_log_q_gain = 0.0L;
    long double dynamic_maximum_local_q = 0.0L;
    long double dynamic_maximum_nonlocal_q = 0.0L;
    long double dynamic_maximum_near_nonlocal_q = 0.0L;
    long double dynamic_maximum_far_nonlocal_q = 0.0L;
    long double dynamic_maximum_selected_gap_tail_q = 0.0L;
    long double dynamic_q_log_growth_ratio = 0.0L;
    long double dynamic_q_derivative_error = 0.0L;
    long double strong_l4_envelope = 0.0L;
    long double envelope_utilization = 0.0L;
    long double dynamic_maximum_enstrophy = 0.0L;
    long double dynamic_maximum_vorticity = 0.0L;
    long double dynamic_maximum_holder_half = 0.0L;
    long double dynamic_maximum_stretch_alignment = 0.0L;
    long double dynamic_nonlocal_vortex_fraction = 0.0L;
    long double dynamic_partition_residual = 0.0L;
    long double dynamic_final_energy = 0.0L;
    long double dynamic_energy_balance_residual = 0.0L;
    long double dynamic_integral_absolute_local_vortex = 0.0L;
    long double dynamic_integral_absolute_nonlocal_vortex = 0.0L;
    long double dynamic_integral_absolute_total_vortex = 0.0L;
    int dynamic_geometry_samples = 0;
    int dynamic_evaluations = 0;
    int dynamic_winning_restart = 0;
    int dynamic_accepted_mutations = 0;
    int dynamic_accepted_gradient_steps = 0;
    long double dynamic_sobolev_value = 0.0L;
    std::vector<long double> dynamic_restart_objectives;
    std::vector<AdversaryGradientTracePoint> dynamic_gradient_trace;
};

struct AdversaryReport {
    int workers = 1;
    std::string backend;
    std::string dynamic_objective;
    std::string dynamic_optimizer;
    std::string gradient_method;
    int minimum_dyadic_gap = 2;
    int sobolev_order = 0;
    long double sobolev_cap = 0.0L;
    long double critical_density_shift = 0.0L;
    int restarts = 0;
    int dynamic_restarts = 1;
    bool dynamic_replay_each_cutoff = false;
    int generations = 0;
    int dynamic_generations = 0;
    long double mutation = 0.0L;
    std::uint64_t seed = 0;
    long double viscosity = 0.0L;
    long double time = 0.0L;
    long double requested_dt = 0.0L;
    long double q_growth_ratio = 1.0L;
    long double q_cutoff_log_slope = 0.0L;
    bool embedding_monotonicity = false;
    std::vector<AdversaryReportRow> rows;
};

class AdversaryReporter {
public:
    static void write_console(const AdversaryReport& report, std::ostream& out);
    static void write_json(const AdversaryReport& report, std::ostream& out);
};

}  // namespace lemma
