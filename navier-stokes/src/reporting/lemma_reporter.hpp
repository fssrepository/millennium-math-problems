#pragma once

#include "dyadic_shell_bounds.hpp"
#include "far_tail_closure.hpp"
#include "helical_triad_ledger.hpp"
#include "periodic_shell_geometry.hpp"
#include "periodic_tail_bound.hpp"
#include "transition_block_scaling.hpp"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>

namespace lemma {

struct LemmaReport {
    std::size_t candidate_count = 0;
    std::string minimum_young_power;
    std::string minimizer_energy;
    std::string minimizer_enstrophy;
    std::string minimizer_palinstrophy;
    std::string young_multiplier_power;
    std::string pointwise_depletion_power;
    std::string energy_depletion_power;
    bool universal_quarter_depletion = false;
    bool closing_candidate_exists = false;
    std::string fixed_energy_q_exponent;
    bool pointwise_q_scale_compatible = false;
    std::string critical_density_exponent;
    std::string time_exponent;
    std::string integrated_l4_exponent;
    bool integrated_l4_scale_critical = false;
    bool exact_strong_l4_factorization = false;
    bool uniform_q_closes_l4 = false;
    std::string dyadic_advecting_gap_decay;
    std::string dyadic_target_gap_decay;
    std::string dyadic_l4_density_gap_decay;
    std::string dyadic_l4_density_enstrophy_power;
    bool dyadic_tail_summable = false;
    bool dyadic_energy_closes_time_integral = false;
    std::string dyadic_post_young_gap_decay;
    std::string dyadic_post_young_enstrophy_power;
    std::string moving_gap_log_enstrophy_slope;
    std::string moving_gap_remaining_enstrophy_power;
    bool moving_gap_closes_far_tail = false;
    DyadicShellRandomCertificate dyadic_shell_bounds;
    PeriodicShellGeometryCertificate periodic_shell_geometry;
    PeriodicTailBoundCertificate periodic_tail_bound;
    FarTailClosureCertificate far_tail_closure;
    TransitionBlockScalingReport transition_block_scaling;
    HelicalTriadCertificate helical_triad_certificate;
    long double helical_adversary_initial_objective = 0.0L;
    long double helical_adversary_final_objective = 0.0L;
    int helical_adversary_accepted_steps = 0;
    int helical_adversary_evaluations = 0;
    int helical_adversary_restarts = 0;
    int helical_adversary_threads = 0;
    long double helical_trajectory_initial_objective = 0.0L;
    long double helical_trajectory_final_objective = 0.0L;
    int helical_trajectory_accepted_steps = 0;
    int helical_trajectory_evaluations = 0;
    int helical_trajectory_restarts = 0;
    int helical_trajectory_threads = 0;
    int triad_cutoff = 0;
    int triad_modes = 0;
    int triad_samples = 0;
    std::uint64_t seed = 0;
    long double energy_residual = 0.0L;
    long double divergence_residual = 0.0L;
    long double reality_residual = 0.0L;
    long double classical_ratio = 0.0L;
    long double detailed_triad_residual = 0.0L;
    long double relative_detailed_triad_residual = 0.0L;
    long double nonlocal_absolute_fraction = 0.0L;
    long double flux_efficiency = 0.0L;
    long double local_cumulative_flux = 0.0L;
    long double nonlocal_cumulative_flux = 0.0L;
    long double flux_partition_residual = 0.0L;
    bool nonzero_vortex_stretching = false;
};

class LemmaReporter {
public:
    static void write_console(const LemmaReport& report, std::ostream& out);
    static void write_json(const LemmaReport& report, std::ostream& out);
};

}  // namespace lemma
