#include "lemma_engine.hpp"
#include "adversary_reporter.hpp"
#include "dyadic_shell_bounds.hpp"
#include "dynamic_adversary.hpp"
#include "family_reporter.hpp"
#include "far_tail_closure.hpp"
#include "gradient_adversary.hpp"
#include "helical_triad_ledger.hpp"
#include "helical_gap_ledger.hpp"
#include "helical_sector_objective.hpp"
#include "helical_sector_adversary.hpp"
#include "helical_sector_adjoint.hpp"
#include "helical_trajectory_adversary.hpp"
#include "proof_scaling.hpp"
#include "parallel_executor.hpp"
#include "periodic_shell_geometry.hpp"
#include "periodic_tail_bound.hpp"
#include "orthogonal_triad_geometry.hpp"
#include "lemma_adversary.hpp"
#include "lemma_reporter.hpp"
#include "local_signature_geometry.hpp"
#include "local_signature_density.hpp"
#include "local_signature_objective.hpp"
#include "local_critical_derivative_ledger.hpp"
#include "local_quartic_identity_ledger.hpp"
#include "local_quartic_commutator.hpp"
#include "local_quartic_closure_target.hpp"
#include "local_quartic_closure_objective.hpp"
#include "local_quartic_projected_residual.hpp"
#include "local_quartic_reduced_ledger.hpp"
#include "local_quartic_shell_ledger.hpp"
#include "local_quartic_shell_envelope.hpp"
#include "local_sld_cyclic_ansatz.hpp"
#include "local_sld_cyclic_krylov_ansatz.hpp"
#include "local_sld_cyclic_trajectory_ansatz.hpp"
#include "local_sld_response_family.hpp"
#include "local_sld_response_hierarchy.hpp"
#include "local_sld_block_objective.hpp"
#include "local_sld_signature_block.hpp"
#include "local_sld_trajectory_adjoint.hpp"
#include "local_sld_trajectory_evaluator.hpp"
#include "local_triad_symmetrizer.hpp"
#include "moving_gap_controller.hpp"
#include "projective_family.hpp"
#include "spectral_adjoint.hpp"
#include "spectral_dynamics.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_objective.hpp"
#include "spectral_state.hpp"
#include "shifted_critical_density.hpp"
#include "shifted_critical_density_budget.hpp"
#include "state_analysis.hpp"
#include "trajectory_analyzer.hpp"
#include "transition_block_scaling.hpp"
#include "triad_commutator.hpp"
#include "triad_tail_envelope.hpp"
#include "triad_verifier.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lemma {
namespace {

using Real = SpectralReal;
using Complex = SpectralComplex;
SpectralGalerkin active_galerkin;
SpectralDynamics active_dynamics(active_galerkin);
SpectralObjective active_objective(active_dynamics);
TrajectoryAnalyzer active_trajectory_analyzer(
    active_galerkin, active_dynamics, active_objective);
SpectralAdjoint active_adjoint(active_dynamics, active_objective);
GradientAdversary active_gradient_adversary(
    active_dynamics, active_objective, active_adjoint);

struct AdversaryResult {
    int cutoff = 0;
    int modes = 0;
    StaticObjective objective;
    SpectralState state;
    int accepted_mutations = 0;
    int evaluations = 0;
};

AdversaryResult optimize_static_depletion(int cutoff, int restarts, int generations,
                                          Real mutation, std::uint64_t seed,
                                          const SpectralState* warm_start = nullptr) {
    if (restarts < 1 || generations < 1) {
        throw std::invalid_argument("adversary restarts and generations must be positive");
    }
    if (!(mutation > 0.0L)) {
        throw std::invalid_argument("--mutation must be positive");
    }
    std::mt19937_64 generator(seed ^
                              (static_cast<std::uint64_t>(cutoff) * 0x9e3779b97f4a7c15ULL));
    AdversaryResult global;
    global.cutoff = cutoff;
    for (int restart = 0; restart < restarts; ++restart) {
        SpectralState current = restart == 0 && warm_start != nullptr
                                    ? SpectralStateFactory::lift(
                                          *warm_start, cutoff, generator)
                                    : SpectralStateFactory::random(cutoff, generator);
        SpectralStateOps::normalize_energy(current);
        StaticObjective current_objective = active_trajectory_analyzer.evaluate_static(current);
        ++global.evaluations;
        if (global.state.waves.empty() ||
            current_objective.energy_level_quantity >
                global.objective.energy_level_quantity) {
            global.state = current;
            global.objective = current_objective;
        }

        for (int generation = 0; generation < generations; ++generation) {
            const Real progress = static_cast<Real>(generation) /
                                  static_cast<Real>(std::max(1, generations - 1));
            const Real scheduled_mutation = mutation * (1.0L - 0.85L * progress);
            SpectralState candidate =
                SpectralStateFactory::mutate(
                    current, scheduled_mutation, generator,
                    generation % 3 != 0);
            const StaticObjective candidate_objective =
                active_trajectory_analyzer.evaluate_static(candidate);
            ++global.evaluations;
            if (candidate_objective.energy_level_quantity >
                current_objective.energy_level_quantity) {
                current = std::move(candidate);
                current_objective = candidate_objective;
                ++global.accepted_mutations;
                if (current_objective.energy_level_quantity >
                    global.objective.energy_level_quantity) {
                    global.state = current;
                    global.objective = current_objective;
                }
            }
        }
    }
    global.modes = static_cast<int>(global.state.waves.size());
    return global;
}

AdversaryResult optimize_static_depletion_parallel(
    int cutoff, int restarts, int generations, Real mutation, std::uint64_t seed,
    const SpectralState* warm_start, const LemmaAdversary& adversary) {
    std::mt19937_64 layout_generator(0);
    const SpectralState layout =
        SpectralStateFactory::random(cutoff, layout_generator);
    static_cast<void>(SpectralStateOps::interactions(layout));

    std::vector<AdversaryResult> partial(static_cast<std::size_t>(restarts));
    adversary.run_restarts(partial.size(), [&](std::size_t restart) {
        const std::uint64_t restart_seed =
            seed + static_cast<std::uint64_t>(restart) * 0x94d049bb133111ebULL;
        partial[restart] = optimize_static_depletion(
            cutoff, 1, generations, mutation, restart_seed,
            restart == 0 ? warm_start : nullptr);
    });

    AdversaryResult result = partial.front();
    int total_evaluations = 0;
    int total_accepted = 0;
    for (const auto& candidate : partial) {
        total_evaluations += candidate.evaluations;
        total_accepted += candidate.accepted_mutations;
        if (candidate.objective.energy_level_quantity >
            result.objective.energy_level_quantity) {
            result = candidate;
        }
    }
    result.evaluations = total_evaluations;
    result.accepted_mutations = total_accepted;
    return result;
}

void write_spectral_state(const std::string& path, const AdversaryResult& result) {
    std::ostringstream metadata;
    metadata << std::setprecision(20)
             << "cutoff=" << result.cutoff
             << " energy=" << static_cast<double>(result.objective.energy)
             << " Q=D^4*Z="
             << static_cast<double>(
                    result.objective.energy_level_quantity);
    SpectralStateWriter::write_tsv(path, result.state, metadata.str());
}

}  // namespace

int run(const Options& options, std::ostream& out) {
    const ScalingCertificate scaling =
        ScalingAnalyzer::analyze_monomials(options.exponent_denominator);
    const ConcentrationScaling concentration =
        ScalingAnalyzer::analyze_concentration();
    const StrongL4Reduction strong_l4 =
        ScalingAnalyzer::analyze_strong_l4_reduction();
    const DyadicTailScaling dyadic_tail =
        ScalingAnalyzer::analyze_dyadic_tail();
    const DyadicShellRandomCertificate dyadic_shell_bounds =
        DyadicShellBounds::verify_random(
            32, 2, 512, options.seed ^ UINT64_C(0xd1a61c5e11));
    const PeriodicShellGeometryCertificate periodic_shell_geometry =
        PeriodicShellGeometry::certify(
            5, 512, options.seed ^ UINT64_C(0x5e1106e7));
    const PeriodicTailBoundCertificate periodic_tail_bound =
        PeriodicTailBound::verify_random(
            3, 2, 4, options.seed ^ UINT64_C(0xfa47a11),
            periodic_shell_geometry);
    const FarTailClosureCertificate far_tail_closure =
        FarTailClosure::verify_random(
            0.1L, 2, 512, options.seed ^ UINT64_C(0xc105ed),
            periodic_shell_geometry);
    const TransitionBlockScalingReport transition_block_scaling =
        TransitionBlockScaling::analyze();
    const HelicalTriadCertificate helical_triad_certificate =
        HelicalTriadLedger::verify_random(
            2, 4, options.seed ^ UINT64_C(0x4e11ca1));
    HelicalSectorAdversaryOptions helical_adversary_options;
    helical_adversary_options.iterations = 8;
    helical_adversary_options.line_search_steps = 16;
    helical_adversary_options.initial_step = 0.1L;
    constexpr int helical_restart_count = 12;
    const LemmaAdversary helical_restart_executor(12);
    std::mt19937_64 helical_layout_generator(0);
    const SpectralState helical_layout =
        SpectralStateFactory::random(2, helical_layout_generator);
    static_cast<void>(SpectralStateOps::interactions(helical_layout));
    std::array<HelicalSectorAdversaryResult, helical_restart_count>
        helical_restart_results;
    helical_restart_executor.run_restarts(
        helical_restart_results.size(), [&](std::size_t restart) {
            std::mt19937_64 generator(
                (options.seed ^ UINT64_C(0xad6e25a)) +
                static_cast<std::uint64_t>(restart) *
                    UINT64_C(0x9e3779b97f4a7c15));
            SpectralState state = SpectralStateFactory::random(2, generator);
            SpectralStateOps::normalize_energy(state);
            helical_restart_results[restart] =
                HelicalSectorAdversary::maximize(
                    state, helical_adversary_options);
        });
    HelicalSectorAdversaryResult helical_adversary =
        helical_restart_results.front();
    int helical_total_evaluations = 0;
    for (const HelicalSectorAdversaryResult& candidate :
         helical_restart_results) {
        helical_total_evaluations += candidate.evaluations;
        if (candidate.objective > helical_adversary.objective) {
            helical_adversary = candidate;
        }
    }
    const HelicalSectorAdjoint helical_trajectory_adjoint(active_dynamics);
    HelicalTrajectoryAdversaryOptions helical_trajectory_options;
    helical_trajectory_options.iterations = 4;
    helical_trajectory_options.line_search_steps = 16;
    helical_trajectory_options.trajectory_steps = 2;
    helical_trajectory_options.initial_step = 0.1L;
    helical_trajectory_options.viscosity = 0.1L;
    helical_trajectory_options.time_step = 0.001L;
    std::array<HelicalTrajectoryAdversaryResult, helical_restart_count>
        helical_trajectory_results;
    helical_restart_executor.run_restarts(
        helical_trajectory_results.size(), [&](std::size_t restart) {
            std::mt19937_64 generator(
                (options.seed ^ UINT64_C(0x7a6ec70)) +
                static_cast<std::uint64_t>(restart) *
                    UINT64_C(0xbf58476d1ce4e5b9));
            SpectralState state = SpectralStateFactory::random(2, generator);
            SpectralStateOps::normalize_energy(state);
            helical_trajectory_results[restart] =
                HelicalTrajectoryAdversary::maximize(
                    state, helical_trajectory_options,
                    helical_trajectory_adjoint);
        });
    HelicalTrajectoryAdversaryResult helical_trajectory =
        helical_trajectory_results.front();
    int helical_trajectory_total_evaluations = 0;
    for (const HelicalTrajectoryAdversaryResult& candidate :
         helical_trajectory_results) {
        helical_trajectory_total_evaluations += candidate.evaluations;
        if (candidate.objective > helical_trajectory.objective) {
            helical_trajectory = candidate;
        }
    }
    const TriadCertificate triads =
        TriadVerifier::analyze(
            options.triad_cutoff, options.triad_samples, options.seed);
    LemmaReport report;
    report.candidate_count = scaling.candidates.size();
    report.minimum_young_power = scaling.minimum_young_power.str();
    report.minimizer_energy = scaling.minimizer.energy.str();
    report.minimizer_enstrophy = scaling.minimizer.enstrophy.str();
    report.minimizer_palinstrophy = scaling.minimizer.palinstrophy.str();
    report.young_multiplier_power =
        scaling.minimizer.young_multiplier_power.str();
    report.pointwise_depletion_power =
        scaling.minimizer.pointwise_linear_depletion_power.str();
    report.energy_depletion_power =
        scaling.minimizer.energy_integrable_depletion_power.str();
    report.universal_quarter_depletion = scaling.universal_quarter_depletion;
    report.closing_candidate_exists = scaling.closing_candidate_exists;
    report.fixed_energy_q_exponent =
        concentration.fixed_energy_pointwise_q.str();
    report.pointwise_q_scale_compatible =
        concentration.pointwise_candidate_scale_compatible;
    report.critical_density_exponent =
        concentration.natural_critical_integrand.str();
    report.time_exponent = concentration.time.str();
    report.integrated_l4_exponent = concentration.natural_integrated_l4.str();
    report.integrated_l4_scale_critical =
        concentration.integrated_candidate_scale_critical;
    report.exact_strong_l4_factorization =
        strong_l4.exact_density_factorization;
    report.uniform_q_closes_l4 =
        strong_l4.closes_integrated_l4_from_uniform_q;
    report.dyadic_advecting_gap_decay =
        dyadic_tail.low_advecting_gap_decay.str();
    report.dyadic_target_gap_decay =
        dyadic_tail.low_target_gap_decay.str();
    report.dyadic_l4_density_gap_decay =
        dyadic_tail.l4_density_gap_decay.str();
    report.dyadic_l4_density_enstrophy_power =
        dyadic_tail.l4_density_enstrophy_power.str();
    report.dyadic_tail_summable =
        dyadic_tail.frequency_tail_is_summable;
    report.dyadic_energy_closes_time_integral =
        dyadic_tail.energy_identity_closes_time_integral;
    report.dyadic_post_young_gap_decay =
        dyadic_tail.post_young_gap_decay.str();
    report.dyadic_post_young_enstrophy_power =
        dyadic_tail.post_young_enstrophy_power.str();
    report.moving_gap_log_enstrophy_slope =
        dyadic_tail.moving_gap_log_enstrophy_slope.str();
    report.moving_gap_remaining_enstrophy_power =
        dyadic_tail.moving_gap_remaining_enstrophy_power.str();
    report.moving_gap_closes_far_tail =
        dyadic_tail.moving_gap_closes_far_tail;
    report.dyadic_shell_bounds = dyadic_shell_bounds;
    report.periodic_shell_geometry = periodic_shell_geometry;
    report.periodic_tail_bound = periodic_tail_bound;
    report.far_tail_closure = far_tail_closure;
    report.transition_block_scaling = transition_block_scaling;
    report.helical_triad_certificate = helical_triad_certificate;
    report.helical_adversary_initial_objective =
        helical_adversary.initial_objective;
    report.helical_adversary_final_objective = helical_adversary.objective;
    report.helical_adversary_accepted_steps =
        helical_adversary.accepted_steps;
    report.helical_adversary_restarts = helical_restart_count;
    report.helical_adversary_threads = helical_restart_executor.threads();
    report.helical_adversary_evaluations = helical_total_evaluations;
    report.helical_trajectory_initial_objective =
        helical_trajectory.initial_objective;
    report.helical_trajectory_final_objective =
        helical_trajectory.objective;
    report.helical_trajectory_accepted_steps =
        helical_trajectory.accepted_steps;
    report.helical_trajectory_evaluations =
        helical_trajectory_total_evaluations;
    report.helical_trajectory_restarts = helical_restart_count;
    report.helical_trajectory_threads = helical_restart_executor.threads();
    report.triad_cutoff = options.triad_cutoff;
    report.triad_modes = triads.modes;
    report.triad_samples = triads.samples;
    report.seed = options.seed;
    report.energy_residual = triads.maximum_normalized_energy_residual;
    report.divergence_residual = triads.maximum_divergence_residual;
    report.reality_residual = triads.maximum_reality_residual;
    report.classical_ratio = triads.maximum_classical_ratio;
    report.detailed_triad_residual = triads.maximum_detailed_triad_residual;
    report.relative_detailed_triad_residual =
        triads.maximum_relative_detailed_triad_residual;
    report.nonlocal_absolute_fraction =
        triads.maximum_nonlocal_absolute_fraction;
    report.flux_efficiency = triads.maximum_flux_efficiency;
    report.local_cumulative_flux = triads.maximum_local_cumulative_flux;
    report.nonlocal_cumulative_flux = triads.maximum_nonlocal_cumulative_flux;
    report.flux_partition_residual = triads.maximum_flux_partition_residual;
    report.nonzero_vortex_stretching = triads.nonzero_vortex_stretching_seen;
    LemmaReporter::write_console(report, out);

    if (!options.certificate_path.empty()) {
        std::ofstream certificate(options.certificate_path);
        if (!certificate) {
            throw std::runtime_error("cannot open certificate: " + options.certificate_path);
        }
        LemmaReporter::write_json(report, certificate);
        out << "Certificate written to " << options.certificate_path << '\n';
    }

    const bool passed = scaling.has_absorbable_candidate &&
                        scaling.minimum_young_power == Rational(3) &&
                        scaling.universal_quarter_depletion &&
                        !concentration.pointwise_candidate_scale_compatible &&
                        concentration.integrated_candidate_scale_critical &&
                        strong_l4.exact_density_factorization &&
                        strong_l4.closes_integrated_l4_from_uniform_q &&
                        dyadic_tail.frequency_tail_is_summable &&
                        !dyadic_tail.energy_identity_closes_time_integral &&
                        dyadic_tail.moving_gap_closes_far_tail &&
                        dyadic_shell_bounds.all_bounds_hold &&
                        periodic_shell_geometry.all_bounds_hold &&
                        periodic_tail_bound.all_bounds_hold &&
                        far_tail_closure.all_bounds_hold &&
                        !transition_block_scaling
                             .energy_identity_closes_transition_block &&
                        helical_triad_certificate
                            .all_reconstruction_checks_hold &&
                        helical_triad_certificate
                            .nonzero_pure_homochiral_local_seen &&
                        helical_adversary.objective >
                            helical_adversary.initial_objective &&
                        helical_trajectory.objective >
                            helical_trajectory.initial_objective &&
                        !scaling.closing_candidate_exists &&
                        triads.maximum_normalized_energy_residual < 1e-15L &&
                        triads.maximum_divergence_residual < 1e-15L &&
                        triads.maximum_reality_residual < 1e-15L &&
                        triads.maximum_detailed_triad_residual < 1e-15L &&
                        triads.maximum_flux_partition_residual < 1e-15L &&
                        triads.nonzero_vortex_stretching_seen;
    return passed ? 0 : 2;
}

bool self_test(std::ostream& out) {
    const ScalingCertificate scaling = ScalingAnalyzer::analyze_monomials(32);
    const ConcentrationScaling concentration =
        ScalingAnalyzer::analyze_concentration();
    const StrongL4Reduction strong_l4 =
        ScalingAnalyzer::analyze_strong_l4_reduction();
    const ShiftedCriticalDensityCertificate shifted_density =
        ShiftedCriticalDensityLemma::analyze();
    const DyadicTailScaling dyadic_tail =
        ScalingAnalyzer::analyze_dyadic_tail();
    const DyadicShellRandomCertificate dyadic_shell_bounds =
        DyadicShellBounds::verify_random(32, 2, 512, 1701);
    const PeriodicShellGeometryCertificate periodic_shell_geometry =
        PeriodicShellGeometry::certify(5, 512, 1702);
    const PeriodicTailBoundCertificate periodic_tail_bound =
        PeriodicTailBound::verify_random(
            3, 2, 4, 1703, periodic_shell_geometry);
    const FarTailClosureCertificate far_tail_closure =
        FarTailClosure::verify_random(
            0.1L, 2, 512, 1704, periodic_shell_geometry);
    const TransitionBlockScalingReport transition_block_scaling =
        TransitionBlockScaling::analyze();
    const bool rational_ok = Rational(1, 2) + Rational(1, 3) == Rational(5, 6) &&
                             Rational(3, 4) * Rational(8, 9) == Rational(2, 3);
    const bool scaling_ok = scaling.has_absorbable_candidate &&
                            scaling.minimum_young_power == Rational(3) &&
                            scaling.universal_quarter_depletion &&
                            !scaling.closing_candidate_exists;
    const bool concentration_ok =
        concentration.fixed_energy_pointwise_q == Rational(2) &&
        !concentration.pointwise_candidate_scale_compatible &&
        concentration.natural_integrated_l4 == Rational(0) &&
        concentration.integrated_candidate_scale_critical;
    const bool strong_l4_ok = strong_l4.exact_density_factorization &&
                              strong_l4.closes_integrated_l4_from_uniform_q;
    const bool shifted_density_ok = shifted_density.shift_matches_density &&
        shifted_density.gronwall_coefficient_is_critical &&
        shifted_density.energy_identity_closes_conditionally &&
        shifted_density.density_amplitude_degree == Rational(4) &&
        shifted_density.density_scaling_exponent == Rational(2);
    const bool dyadic_tail_scaling_ok =
        dyadic_tail.low_advecting_gap_decay == Rational(1, 2) &&
        dyadic_tail.low_advected_gap_decay == Rational(1, 2) &&
        dyadic_tail.low_target_gap_decay == Rational(3, 2) &&
        dyadic_tail.l4_density_gap_decay == Rational(2) &&
        dyadic_tail.l4_density_enstrophy_power == Rational(2) &&
        dyadic_tail.l4_density_palinstrophy_power == Rational(0) &&
        dyadic_tail.young_palinstrophy_conjugate == Rational(4, 3) &&
        dyadic_tail.young_remainder_conjugate == Rational(4) &&
        dyadic_tail.post_young_gap_decay == Rational(2) &&
        dyadic_tail.post_young_enstrophy_power == Rational(3) &&
        dyadic_tail.post_young_inverse_viscosity_power == Rational(3) &&
        dyadic_tail.moving_gap_remaining_enstrophy_power == Rational(1) &&
        dyadic_tail.frequency_tail_is_summable &&
        !dyadic_tail.energy_identity_closes_time_integral &&
        dyadic_tail.moving_gap_closes_far_tail;
    const bool dyadic_shell_bounds_ok =
        dyadic_shell_bounds.all_bounds_hold &&
        dyadic_shell_bounds.maximum_high_moment_ratio <= 1.0L &&
        dyadic_shell_bounds.maximum_low_one_derivative_ratio <= 1.0L &&
        dyadic_shell_bounds.maximum_low_three_derivative_ratio <= 1.0L &&
        dyadic_shell_bounds.maximum_one_gain_tail_ratio <= 1.0L &&
        dyadic_shell_bounds.maximum_three_gain_tail_ratio <= 1.0L;
    const bool periodic_shell_geometry_ok =
        periodic_shell_geometry.all_bounds_hold &&
        periodic_shell_geometry.lattice_count_constant == 64.0L &&
        periodic_shell_geometry.l2_to_linf_bernstein_constant == 8.0L &&
        periodic_shell_geometry.gradient_bernstein_constant == 16.0L &&
        periodic_shell_geometry.separated_high_shell_neighbor_width == 1 &&
        periodic_shell_geometry.maximum_count_ratio <= 1.0L &&
        periodic_shell_geometry.maximum_one_gain_overlap_ratio <= 1.0L &&
        periodic_shell_geometry.maximum_three_gain_overlap_ratio <= 1.0L;
    const bool periodic_tail_bound_ok =
        periodic_tail_bound.all_bounds_hold &&
        periodic_tail_bound.nonzero_tail_seen &&
        periodic_tail_bound.maximum_bound_ratio <= 1.0L;
    const bool far_tail_closure_ok =
        far_tail_closure.all_bounds_hold &&
        far_tail_closure.maximum_normalized_remainder_ratio <= 1.0L;
    const bool transition_block_scaling_ok =
        transition_block_scaling.post_young_logarithm_power == Rational(4) &&
        transition_block_scaling.post_young_enstrophy_power == Rational(3) &&
        transition_block_scaling.required_pointwise_depletion_power ==
            Rational(1, 2) &&
        !transition_block_scaling
             .logarithmic_band_count_changes_polynomial_power &&
        !transition_block_scaling.energy_identity_closes_transition_block &&
        TransitionBlockScaling::normalized_remainder(2.0L) <
            TransitionBlockScaling::normalized_remainder(16.0L) &&
        TransitionBlockScaling::normalized_remainder(16.0L) <
            TransitionBlockScaling::normalized_remainder(256.0L);
    const std::array<Real, 7> moving_gap_enstrophies{
        0.0L, 0.25L, 1.0L, 1.1L, 4.0L, 4.1L, 1024.0L};
    bool moving_gap_controller_ok = true;
    for (const Real enstrophy : moving_gap_enstrophies) {
        const MovingGapDecision decision =
            MovingGapController::decide(enstrophy, 2);
        moving_gap_controller_ok = moving_gap_controller_ok &&
            decision.minimum_gap ==
                2 + decision.logarithmic_gap &&
            decision.normalized_cubic_remainder_ratio <=
                1.0L + 1e-15L &&
            decision.base_weighted_remainder_ratio <=
                0.0625L + 1e-15L;
    }
    moving_gap_controller_ok = moving_gap_controller_ok &&
        MovingGapController::decide(1.1L, 2).minimum_gap == 3 &&
        MovingGapController::decide(4.0L, 2).minimum_gap == 4 &&
        MovingGapController::decide(4.1L, 2).minimum_gap == 5;
    const TriadCertificate triads = TriadVerifier::analyze(2, 2, 7);
    const bool triad_ok = triads.maximum_normalized_energy_residual < 1e-15L &&
                          triads.maximum_divergence_residual < 1e-15L &&
                          triads.maximum_detailed_triad_residual < 1e-15L &&
                          triads.maximum_flux_partition_residual < 1e-15L &&
                          triads.nonzero_vortex_stretching_seen;
    std::mt19937_64 helical_generator(11);
    SpectralState helical_state =
        SpectralStateFactory::random(2, helical_generator);
    SpectralStateOps::normalize_energy(helical_state);
    const HelicalTriadReport helical =
        HelicalTriadLedger::analyze(helical_state);
    const HelicalGapLedgerReport helical_gaps =
        HelicalGapLedger::analyze(helical_state);
    const LocalTriadSymmetryReport local_symmetry =
        LocalTriadSymmetrizer::analyze(helical_state);
    const OrthogonalTriadGeometryCertificate orthogonal_geometry =
        OrthogonalTriadGeometry::certify(5);
    const OrthogonalTriadClosure orthogonal_closure =
        OrthogonalTriadGeometry::analyze_closure();
    const LocalSignatureGeometryCertificate local_signature_geometry =
        LocalSignatureGeometry::certify(3);
    const SignatureFamilyClosure signature_family_closure =
        LocalSignatureGeometry::analyze_closure();
    SpectralState positive_helical_state =
        HelicalTriadLedger::project_helicity(helical_state, 1);
    SpectralStateOps::normalize_energy(positive_helical_state);
    const HelicalTriadReport positive_helical =
        HelicalTriadLedger::analyze(positive_helical_state);
    SpectralState negative_helical_state =
        HelicalTriadLedger::project_helicity(helical_state, -1);
    SpectralStateOps::normalize_energy(negative_helical_state);
    const HelicalTriadReport negative_helical =
        HelicalTriadLedger::analyze(negative_helical_state);
    const bool helical_ok =
        helical.relative_velocity_reconstruction_residual < 1e-15L &&
        helical.relative_total_reconstruction_residual < 1e-15L &&
        helical.relative_local_reconstruction_residual < 1e-15L &&
        std::abs(helical.signed_local_stretching -
                 helical.homochiral_local_stretching -
                 helical.heterochiral_local_stretching) < 1e-15L;
    const bool helical_gap_ok =
        helical_gaps.maximum_gap_reconstruction_residual < 1e-15L &&
        helical_gaps.total_reconstruction_residual < 1e-15L &&
        !helical_gaps.gaps.empty() &&
        std::abs(helical_gaps.gaps.front().homochiral_signed -
                 helical.homochiral_local_stretching) < 1e-15L &&
        std::abs(helical_gaps.gaps.front().heterochiral_signed -
                 helical.heterochiral_local_stretching) < 1e-15L;
    const bool local_symmetry_ok =
        local_symmetry.maximum_energy_cancellation_residual < 1e-15L &&
        local_symmetry.local_reconstruction_residual < 1e-15L &&
        local_symmetry.maximum_frequency_spread_bound_ratio <=
            1.0L + 1e-15L &&
        local_symmetry.coherent_signature_count > 0 &&
        local_symmetry.effective_coherent_signature_count >= 1.0L &&
        local_symmetry.effective_coherent_signature_count <=
            static_cast<Real>(local_symmetry.signatures.size()) +
                1e-12L &&
        local_symmetry.dominant_coherent_signature_fraction > 0.0L &&
        local_symmetry.dominant_coherent_signature_fraction <= 1.0L &&
        local_symmetry.signed_signature_cancellation_ratio <=
            1.0L + 1e-15L &&
        local_symmetry.signed_signature_amplification <=
            std::sqrt(local_symmetry.effective_coherent_signature_count) +
                1e-15L &&
        local_symmetry.local_triads > 0;
    const bool orthogonal_geometry_ok =
        orthogonal_geometry.all_degree_bounds_hold &&
        orthogonal_geometry.maximum_input_degree_ratio <= 1.0L &&
        orthogonal_geometry.maximum_target_degree_ratio <= 1.0L &&
        orthogonal_closure.transfer_frequency_power == Rational(7, 2) &&
        orthogonal_closure.generic_local_transfer_frequency_power ==
            Rational(9, 2) &&
        orthogonal_closure.critical_transfer_frequency_power == Rational(4) &&
        orthogonal_closure.transfer_to_viscosity_frequency_power ==
            Rational(-1, 2) &&
        orthogonal_closure.high_frequency_absorbable_from_energy &&
        orthogonal_closure.orthogonal_degree_is_subcritical &&
        orthogonal_closure.generic_local_degree_is_supercritical;
    const bool local_signature_geometry_ok =
        local_signature_geometry.all_fixed_signature_degree_bounds_hold &&
        local_signature_geometry.maximum_input_degree_ratio <= 1.0L &&
        local_signature_geometry.maximum_target_degree_ratio <= 1.0L &&
        signature_family_closure.finite_signature_family
                .transfer_frequency_power == Rational(7, 2) &&
        signature_family_closure.finite_signature_family
                .energy_level_high_frequency_absorption &&
        signature_family_closure.critical_signature_family
                .transfer_frequency_power == Rational(4) &&
        !signature_family_closure.critical_signature_family
                .energy_level_high_frequency_absorption &&
        signature_family_closure.dense_signature_family
                .transfer_frequency_power == Rational(9, 2) &&
        !signature_family_closure.dense_signature_family
                .energy_level_high_frequency_absorption &&
        signature_family_closure.square_summed_fixed_signature_bound &&
        signature_family_closure.effective_count_replaces_raw_count &&
        signature_family_closure.critical_signed_amplification_power ==
            Rational(1, 2) &&
        signature_family_closure
            .signed_amplification_preserves_cancellation &&
        signature_family_closure
                .closing_requires_sublinear_signature_count;
    const bool pure_helical_ok =
        positive_helical.negative_helical_energy < 1e-15L &&
        negative_helical.positive_helical_energy < 1e-15L &&
        positive_helical.heterochiral_absolute_local_stretching < 1e-15L &&
        negative_helical.heterochiral_absolute_local_stretching < 1e-15L &&
        positive_helical.relative_local_reconstruction_residual < 1e-15L &&
        negative_helical.relative_local_reconstruction_residual < 1e-15L;
    std::mt19937_64 fft_generator(19);
    SpectralState fft_state = SpectralStateFactory::random(2, fft_generator);
    SpectralStateOps::normalize_energy(fft_state);
    const auto direct_advection = active_dynamics.advection_direct(fft_state);
    const auto fft_advection = active_dynamics.advection_fft(fft_state);
    Real fft_error2 = 0.0L;
    Real fft_reference2 = 0.0L;
    for (std::size_t mode = 0; mode < direct_advection.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            fft_error2 += std::norm(direct_advection[mode][component] -
                                    fft_advection[mode][component]);
            fft_reference2 += std::norm(direct_advection[mode][component]);
        }
    }
    const Real fft_relative_error =
        std::sqrt(fft_error2 / std::max(1e-30L, fft_reference2));
    const bool fft_ok = fft_relative_error < 1e-14L;

    std::mt19937_64 adjoint_generator(23);
    SpectralState adjoint_state =
        SpectralStateFactory::random(1, adjoint_generator);
    SpectralState tangent_state =
        SpectralStateFactory::random(1, adjoint_generator);
    SpectralState cotangent_state =
        SpectralStateFactory::random(1, adjoint_generator);
    SpectralStateOps::normalize_energy(adjoint_state);
    SpectralStateOps::normalize_energy(tangent_state);
    SpectralStateOps::normalize_energy(cotangent_state);
    const SpectralIncrement& tangent = tangent_state.velocity;
    const SpectralIncrement& cotangent = cotangent_state.velocity;
    auto increment_inner_product = [](const SpectralIncrement& left,
                                      const SpectralIncrement& right) {
        Real result = 0.0L;
        for (std::size_t mode = 0; mode < left.size(); ++mode) {
            result += std::real(dot_hermitian(left[mode], right[mode]));
        }
        return result;
    };
    auto increment_relative_error = [&](const SpectralIncrement& computed,
                                        const SpectralIncrement& reference) {
        Real error2 = 0.0L;
        Real reference2 = 0.0L;
        for (std::size_t mode = 0; mode < computed.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                error2 +=
                    std::norm(computed[mode][component] -
                              reference[mode][component]);
                reference2 += std::norm(reference[mode][component]);
            }
        }
        return std::sqrt(error2 / std::max(1e-30L, reference2));
    };
    SpectralState helical_direction_state =
        SpectralStateFactory::random(2, helical_generator);
    SpectralStateOps::normalize_energy(helical_direction_state);
    constexpr Real helical_gradient_step = 1e-6L;
    const SpectralState helical_plus = active_dynamics.add_increment(
        helical_state, helical_direction_state.velocity,
        helical_gradient_step);
    const SpectralState helical_minus = active_dynamics.add_increment(
        helical_state, helical_direction_state.velocity,
        -helical_gradient_step);
    const LocalSignatureObjectiveValue local_signature_objective =
        LocalSignatureObjective::evaluate(helical_state);
    const LocalSignatureDensitySample local_signature_density =
        LocalSignatureDensity::evaluate(helical_state);
    const SpectralIncrement local_signature_gradient =
        LocalSignatureObjective::signed_amplification_gradient(
            helical_state);
    const Real local_signature_directional = increment_inner_product(
        local_signature_gradient, helical_direction_state.velocity);
    const Real local_signature_finite_difference =
        (LocalSignatureObjective::evaluate(helical_plus)
             .signed_amplification -
         LocalSignatureObjective::evaluate(helical_minus)
             .signed_amplification) /
        (2.0L * helical_gradient_step);
    const Real local_signature_gradient_error = std::abs(
        local_signature_directional -
        local_signature_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(local_signature_directional),
            std::abs(local_signature_finite_difference)));
    const Real local_signature_objective_error = std::abs(
        local_signature_objective.signed_amplification -
        local_symmetry.signed_signature_amplification) /
        std::max(1e-30L,
                 local_symmetry.signed_signature_amplification);
    const SpectralIncrement local_signature_transfer_gradient =
        LocalSignatureObjective::absolute_signed_transfer_gradient(
            helical_state);
    const Real local_signature_transfer_directional =
        increment_inner_product(
            local_signature_transfer_gradient,
            helical_direction_state.velocity);
    const Real local_signature_transfer_finite_difference =
        (std::abs(LocalSignatureObjective::evaluate(helical_plus)
                      .signed_local_transfer) -
         std::abs(LocalSignatureObjective::evaluate(helical_minus)
                      .signed_local_transfer)) /
        (2.0L * helical_gradient_step);
    const Real local_signature_transfer_gradient_error = std::abs(
        local_signature_transfer_directional -
        local_signature_transfer_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(local_signature_transfer_directional),
            std::abs(local_signature_transfer_finite_difference)));
    const bool local_signature_objective_ok =
        local_signature_objective_error < 1e-15L &&
        local_signature_density.factorization_residual < 1e-15L &&
        local_signature_gradient_error < 1e-9L &&
        local_signature_transfer_gradient_error < 1e-9L;
    const HelicalSectorSelection homochiral_selection =
        HelicalSectorSelection::homochiral();
    const HelicalSectorSelection heterochiral_selection =
        HelicalSectorSelection::heterochiral();
    const HelicalSectorSelection broad_heterochiral_selection =
        heterochiral_selection.with_spread(HelicalLocalSpread::broad);
    const HelicalSectorObjectiveValue homochiral_value =
        HelicalSectorObjective::evaluate(
            helical_state, homochiral_selection);
    const HelicalSectorObjectiveValue heterochiral_value =
        HelicalSectorObjective::evaluate(
            helical_state, heterochiral_selection);
    const SpectralIncrement helical_signed_gradient =
        HelicalSectorObjective::signed_stretching_gradient(
            helical_state, heterochiral_selection);
    const SpectralIncrement helical_critical_gradient =
        HelicalSectorObjective::critical_integrand_gradient(
            helical_state, heterochiral_selection);
    const Real helical_signed_directional = increment_inner_product(
        helical_signed_gradient, helical_direction_state.velocity);
    const Real helical_signed_finite_difference =
        (HelicalSectorObjective::evaluate(
             helical_plus, heterochiral_selection)
             .signed_local_stretching -
         HelicalSectorObjective::evaluate(
             helical_minus, heterochiral_selection)
             .signed_local_stretching) /
        (2.0L * helical_gradient_step);
    const Real helical_signed_gradient_error = std::abs(
        helical_signed_directional - helical_signed_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(helical_signed_directional),
            std::abs(helical_signed_finite_difference)));
    const Real helical_critical_directional = increment_inner_product(
        helical_critical_gradient, helical_direction_state.velocity);
    const Real helical_critical_finite_difference =
        (HelicalSectorObjective::evaluate(
             helical_plus, heterochiral_selection)
             .critical_integrand -
         HelicalSectorObjective::evaluate(
             helical_minus, heterochiral_selection)
             .critical_integrand) /
        (2.0L * helical_gradient_step);
    const Real helical_critical_gradient_error = std::abs(
        helical_critical_directional - helical_critical_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(helical_critical_directional),
            std::abs(helical_critical_finite_difference)));
    const SpectralIncrement broad_helical_signed_gradient =
        HelicalSectorObjective::signed_stretching_gradient(
            helical_state, broad_heterochiral_selection);
    const Real broad_helical_signed_directional = increment_inner_product(
        broad_helical_signed_gradient, helical_direction_state.velocity);
    const Real broad_helical_signed_finite_difference =
        (HelicalSectorObjective::evaluate(
             helical_plus, broad_heterochiral_selection)
             .signed_local_stretching -
         HelicalSectorObjective::evaluate(
             helical_minus, broad_heterochiral_selection)
             .signed_local_stretching) /
        (2.0L * helical_gradient_step);
    const Real broad_helical_signed_gradient_error = std::abs(
        broad_helical_signed_directional -
        broad_helical_signed_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(broad_helical_signed_directional),
            std::abs(broad_helical_signed_finite_difference)));
    const SpectralIncrement broad_helical_critical_gradient =
        HelicalSectorObjective::critical_integrand_gradient(
            helical_state, broad_heterochiral_selection);
    const Real broad_helical_critical_directional = increment_inner_product(
        broad_helical_critical_gradient, helical_direction_state.velocity);
    const Real broad_helical_critical_finite_difference =
        (HelicalSectorObjective::evaluate(
             helical_plus, broad_heterochiral_selection)
             .critical_integrand -
         HelicalSectorObjective::evaluate(
             helical_minus, broad_heterochiral_selection)
             .critical_integrand) /
        (2.0L * helical_gradient_step);
    const Real broad_helical_critical_gradient_error = std::abs(
        broad_helical_critical_directional -
        broad_helical_critical_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(broad_helical_critical_directional),
            std::abs(broad_helical_critical_finite_difference)));
    const Real helical_sector_partition_error = std::abs(
        homochiral_value.signed_local_stretching +
        heterochiral_value.signed_local_stretching -
        helical.signed_local_stretching) /
        std::max(1e-30L, helical.homochiral_absolute_local_stretching +
                            helical.heterochiral_absolute_local_stretching);
    const bool helical_sector_objective_ok =
        helical_sector_partition_error < 1e-15L &&
        helical_signed_gradient_error < 1e-9L &&
        helical_critical_gradient_error < 1e-9L &&
        broad_helical_signed_gradient_error < 1e-9L &&
        broad_helical_critical_gradient_error < 1e-9L;
    HelicalSectorAdversaryOptions helical_adversary_options;
    helical_adversary_options.iterations = 3;
    helical_adversary_options.line_search_steps = 16;
    helical_adversary_options.initial_step = 0.1L;
    const HelicalSectorAdversaryResult helical_adversary =
        HelicalSectorAdversary::maximize(
            helical_state, helical_adversary_options);
    const bool helical_adversary_ok =
        helical_adversary.accepted_steps > 0 &&
        helical_adversary.objective > helical_adversary.initial_objective &&
        std::abs(SpectralStateOps::energy(helical_adversary.state) -
                 SpectralStateOps::energy(helical_state)) < 1e-15L;
    constexpr Real helical_trajectory_viscosity = 0.1L;
    constexpr Real helical_trajectory_dt = 0.001L;
    constexpr int helical_trajectory_steps = 2;
    const HelicalSectorAdjoint helical_sector_adjoint(active_dynamics);
    const HelicalSectorTrajectoryGradient helical_trajectory_gradient =
        helical_sector_adjoint.critical_integral_gradient(
            helical_state, helical_trajectory_viscosity,
            helical_trajectory_dt, helical_trajectory_steps,
            heterochiral_selection);
    auto helical_trajectory_integral = [&](SpectralState state) {
        Real integral = 0.5L * helical_trajectory_dt *
            HelicalSectorObjective::evaluate(
                state, heterochiral_selection).critical_integrand;
        for (int step = 0; step < helical_trajectory_steps; ++step) {
            active_dynamics.rk4_step(
                state, helical_trajectory_viscosity,
                helical_trajectory_dt);
            const Real weight = step + 1 == helical_trajectory_steps
                ? 0.5L * helical_trajectory_dt
                : helical_trajectory_dt;
            integral += weight * HelicalSectorObjective::evaluate(
                state, heterochiral_selection).critical_integrand;
        }
        return integral;
    };
    const Real helical_trajectory_directional = increment_inner_product(
        helical_trajectory_gradient.initial_gradient,
        helical_direction_state.velocity);
    const Real helical_trajectory_finite_difference =
        (helical_trajectory_integral(helical_plus) -
         helical_trajectory_integral(helical_minus)) /
        (2.0L * helical_gradient_step);
    const Real helical_trajectory_gradient_error = std::abs(
        helical_trajectory_directional -
        helical_trajectory_finite_difference) /
        std::max(1e-30L, std::max(
            std::abs(helical_trajectory_directional),
            std::abs(helical_trajectory_finite_difference)));
    const bool helical_trajectory_adjoint_ok =
        helical_trajectory_gradient_error < 1e-9L;
    HelicalTrajectoryAdversaryOptions helical_trajectory_options;
    helical_trajectory_options.iterations = 2;
    helical_trajectory_options.line_search_steps = 16;
    helical_trajectory_options.trajectory_steps = helical_trajectory_steps;
    helical_trajectory_options.initial_step = 0.1L;
    helical_trajectory_options.viscosity = helical_trajectory_viscosity;
    helical_trajectory_options.time_step = helical_trajectory_dt;
    const HelicalTrajectoryAdversaryResult helical_trajectory_adversary =
        HelicalTrajectoryAdversary::maximize(
            helical_state, helical_trajectory_options,
            helical_sector_adjoint);
    const bool helical_trajectory_adversary_ok =
        helical_trajectory_adversary.accepted_steps > 0 &&
        helical_trajectory_adversary.objective >
            helical_trajectory_adversary.initial_objective &&
        std::abs(SpectralStateOps::energy(
                     helical_trajectory_adversary.state) -
                 SpectralStateOps::energy(helical_state)) < 1e-15L;
    SpectralState fft_tangent_state =
        SpectralStateFactory::random(2, adjoint_generator);
    SpectralState fft_cotangent_state =
        SpectralStateFactory::random(2, adjoint_generator);
    SpectralStateOps::normalize_energy(fft_tangent_state);
    SpectralStateOps::normalize_energy(fft_cotangent_state);
    const SpectralIncrement fft_jvp_direct =
        active_dynamics.advection_jvp_direct(
            fft_state, fft_tangent_state.velocity);
    const SpectralIncrement fft_jvp = active_dynamics.advection_jvp_fft(
        fft_state, fft_tangent_state.velocity);
    const SpectralIncrement fft_vjp_direct =
        active_dynamics.advection_vjp_direct(
            fft_state, fft_cotangent_state.velocity);
    const SpectralIncrement fft_vjp = active_dynamics.advection_vjp_fft(
        fft_state, fft_cotangent_state.velocity);
    const Real fft_jvp_oracle_error =
        increment_relative_error(fft_jvp, fft_jvp_direct);
    const Real fft_vjp_oracle_error =
        increment_relative_error(fft_vjp, fft_vjp_direct);
    const Real fft_duality_left = increment_inner_product(
        fft_cotangent_state.velocity, fft_jvp);
    const Real fft_duality_right = increment_inner_product(
        fft_vjp, fft_tangent_state.velocity);
    const Real fft_adjoint_duality_error =
        std::abs(fft_duality_left - fft_duality_right) /
        std::max(1e-30L,
                 std::max(std::abs(fft_duality_left),
                          std::abs(fft_duality_right)));
    const bool fft_adjoint_ok = fft_jvp_oracle_error < 1e-14L &&
                                fft_vjp_oracle_error < 1e-14L &&
                                fft_adjoint_duality_error < 1e-12L;
    const Real adjoint_viscosity = 0.1L;
    const Real adjoint_dt = 0.001L;
    const SpectralIncrement rhs_tangent = active_dynamics.rhs_jvp(
        adjoint_state, tangent, adjoint_viscosity);
    const SpectralIncrement rhs_cotangent = active_dynamics.rhs_vjp(
        adjoint_state, cotangent, adjoint_viscosity);
    const Real rhs_duality_left =
        increment_inner_product(cotangent, rhs_tangent);
    const Real rhs_duality_right =
        increment_inner_product(rhs_cotangent, tangent);
    const Real rhs_duality_error =
        std::abs(rhs_duality_left - rhs_duality_right) /
        std::max(1e-30L,
                 std::max(std::abs(rhs_duality_left),
                          std::abs(rhs_duality_right)));

    const Real finite_difference_step = 1e-6L;
    const SpectralState rhs_plus_state = active_dynamics.add_increment(
        adjoint_state, tangent, finite_difference_step);
    const SpectralState rhs_minus_state = active_dynamics.add_increment(
        adjoint_state, tangent, -finite_difference_step);
    const SpectralIncrement rhs_plus =
        active_dynamics.rhs(rhs_plus_state, adjoint_viscosity);
    const SpectralIncrement rhs_minus =
        active_dynamics.rhs(rhs_minus_state, adjoint_viscosity);
    SpectralIncrement rhs_finite_difference = rhs_plus;
    for (std::size_t mode = 0; mode < rhs_finite_difference.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            rhs_finite_difference[mode][component] =
                (rhs_plus[mode][component] - rhs_minus[mode][component]) /
                (2.0L * finite_difference_step);
        }
    }
    const Real rhs_jvp_error =
        increment_relative_error(rhs_tangent, rhs_finite_difference);

    const SpectralIncrement rk4_tangent = active_dynamics.rk4_jvp(
        adjoint_state, tangent, adjoint_viscosity, adjoint_dt);
    const SpectralIncrement rk4_cotangent = active_dynamics.rk4_vjp(
        adjoint_state, cotangent, adjoint_viscosity, adjoint_dt);
    const Real rk4_duality_left =
        increment_inner_product(cotangent, rk4_tangent);
    const Real rk4_duality_right =
        increment_inner_product(rk4_cotangent, tangent);
    const Real rk4_duality_error =
        std::abs(rk4_duality_left - rk4_duality_right) /
        std::max(1e-30L,
                 std::max(std::abs(rk4_duality_left),
                          std::abs(rk4_duality_right)));
    SpectralState rk4_plus = rhs_plus_state;
    SpectralState rk4_minus = rhs_minus_state;
    active_dynamics.rk4_step(rk4_plus, adjoint_viscosity, adjoint_dt);
    active_dynamics.rk4_step(rk4_minus, adjoint_viscosity, adjoint_dt);
    SpectralIncrement rk4_finite_difference = rk4_plus.velocity;
    for (std::size_t mode = 0; mode < rk4_finite_difference.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            rk4_finite_difference[mode][component] =
                (rk4_plus.velocity[mode][component] -
                 rk4_minus.velocity[mode][component]) /
                (2.0L * finite_difference_step);
        }
    }
    const Real rk4_jvp_error =
        increment_relative_error(rk4_tangent, rk4_finite_difference);
    const bool adjoint_ok = rhs_jvp_error < 1e-10L &&
                            rhs_duality_error < 1e-12L &&
                            rk4_jvp_error < 1e-10L &&
                            rk4_duality_error < 1e-12L;
    const SpectralIncrement q_gradient =
        active_objective.energy_level_gradient(adjoint_state);
    const Real q_directional_adjoint =
        increment_inner_product(q_gradient, tangent);
    const Real q_plus =
        active_objective.evaluate(rhs_plus_state).energy_level_quantity;
    const Real q_minus =
        active_objective.evaluate(rhs_minus_state).energy_level_quantity;
    const Real q_directional_finite_difference =
        (q_plus - q_minus) / (2.0L * finite_difference_step);
    const Real q_gradient_error =
        std::abs(q_directional_adjoint -
                 q_directional_finite_difference) /
        std::max(1e-30L,
                 std::max(std::abs(q_directional_adjoint),
                          std::abs(q_directional_finite_difference)));
    const bool q_gradient_ok = q_gradient_error < 1e-10L;
    constexpr int trajectory_steps = 3;
    const QTrajectoryGradient trajectory_gradient =
        active_adjoint.terminal_q_gradient(
            adjoint_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps);
    const Real trajectory_directional_adjoint = increment_inner_product(
        trajectory_gradient.initial_gradient, tangent);
    SpectralState trajectory_plus = rhs_plus_state;
    SpectralState trajectory_minus = rhs_minus_state;
    for (int step = 0; step < trajectory_steps; ++step) {
        active_dynamics.rk4_step(
            trajectory_plus, adjoint_viscosity, adjoint_dt);
        active_dynamics.rk4_step(
            trajectory_minus, adjoint_viscosity, adjoint_dt);
    }
    const Real trajectory_q_plus = active_objective
        .evaluate(trajectory_plus)
        .energy_level_quantity;
    const Real trajectory_q_minus = active_objective
        .evaluate(trajectory_minus)
        .energy_level_quantity;
    const Real trajectory_directional_finite_difference =
        (trajectory_q_plus - trajectory_q_minus) /
        (2.0L * finite_difference_step);
    const Real trajectory_gradient_error =
        std::abs(trajectory_directional_adjoint -
                 trajectory_directional_finite_difference) /
        std::max(1e-30L,
                 std::max(std::abs(trajectory_directional_adjoint),
                          std::abs(
                              trajectory_directional_finite_difference)));
    const bool trajectory_gradient_ok =
        trajectory_gradient_error < 1e-10L &&
        trajectory_gradient.objective_step == trajectory_steps &&
        trajectory_gradient.checkpoint_count ==
            static_cast<std::size_t>(trajectory_steps + 1);
    const QTrajectoryGradient q_gain_gradient =
        active_adjoint.q_gain_gradient(
            adjoint_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps);
    const Real q_gain_directional_adjoint = increment_inner_product(
        q_gain_gradient.initial_gradient, tangent);
    const Real q_gain_plus = std::log(trajectory_q_plus / q_plus);
    const Real q_gain_minus = std::log(trajectory_q_minus / q_minus);
    const Real q_gain_directional_finite_difference =
        (q_gain_plus - q_gain_minus) /
        (2.0L * finite_difference_step);
    const Real q_gain_gradient_error =
        std::abs(q_gain_directional_adjoint -
                 q_gain_directional_finite_difference) /
        std::max(1e-30L,
                 std::max(std::abs(q_gain_directional_adjoint),
                          std::abs(q_gain_directional_finite_difference)));
    const bool q_gain_gradient_ok = q_gain_gradient_error < 1e-9L;
    const QTrajectoryGradient q_increase_gradient =
        active_adjoint.q_increase_gradient(
            adjoint_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps);
    const Real q_increase_directional_adjoint = increment_inner_product(
        q_increase_gradient.initial_gradient, tangent);
    const Real q_increase_plus = trajectory_q_plus - q_plus;
    const Real q_increase_minus = trajectory_q_minus - q_minus;
    const Real q_increase_directional_finite_difference =
        (q_increase_plus - q_increase_minus) /
        (2.0L * finite_difference_step);
    const Real q_increase_gradient_error =
        std::abs(q_increase_directional_adjoint -
                 q_increase_directional_finite_difference) /
        std::max(1e-30L,
                 std::max(std::abs(q_increase_directional_adjoint),
                          std::abs(
                              q_increase_directional_finite_difference)));
    const bool q_increase_gradient_ok = q_increase_gradient_error < 1e-9L;
    Real q_increase_divergence_residual = 0.0L;
    Real q_increase_reality_residual = 0.0L;
    for (std::size_t mode = 0;
         mode < q_increase_gradient.initial_gradient.size(); ++mode) {
        const WaveVector wave = adjoint_state.waves[mode];
        const ComplexVector& value =
            q_increase_gradient.initial_gradient[mode];
        q_increase_divergence_residual = std::max(
            q_increase_divergence_residual,
            std::abs(wave_dot(wave, value)));
        const std::size_t negative = adjoint_state.index.at(-wave);
        for (std::size_t component = 0; component < 3; ++component) {
            q_increase_reality_residual = std::max(
                q_increase_reality_residual,
                std::abs(q_increase_gradient.initial_gradient[negative]
                             [component] -
                         std::conj(value[component])));
        }
    }
    const bool q_increase_constraints_ok =
        q_increase_divergence_residual < 1e-15L &&
        q_increase_reality_residual < 1e-15L;
    const QTrajectoryGradient critical_integral_gradient =
        active_adjoint.critical_integral_gradient(
            adjoint_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps);
    const Real critical_integral_directional_adjoint =
        increment_inner_product(
            critical_integral_gradient.initial_gradient, tangent);
    auto discrete_critical_integral = [&](SpectralState state) {
        StaticObjective previous = active_objective.evaluate(state);
        Real integral = 0.0L;
        for (int step = 0; step < trajectory_steps; ++step) {
            active_dynamics.rk4_step(
                state, adjoint_viscosity, adjoint_dt);
            const StaticObjective current = active_objective.evaluate(state);
            integral += 0.5L * adjoint_dt *
                        (previous.critical_integrand +
                         current.critical_integrand);
            previous = current;
        }
        return integral;
    };
    const Real critical_integral_plus =
        discrete_critical_integral(rhs_plus_state);
    const Real critical_integral_minus =
        discrete_critical_integral(rhs_minus_state);
    const Real critical_integral_directional_finite_difference =
        (critical_integral_plus - critical_integral_minus) /
        (2.0L * finite_difference_step);
    const Real critical_integral_gradient_error =
        std::abs(critical_integral_directional_adjoint -
                 critical_integral_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(std::abs(critical_integral_directional_adjoint),
                     std::abs(
                         critical_integral_directional_finite_difference)));
    const bool critical_integral_gradient_ok =
        critical_integral_gradient_error < 1e-9L;
    SpectralState partition_state =
        SpectralStateFactory::random(2, adjoint_generator);
    SpectralState partition_tangent_state =
        SpectralStateFactory::random(2, adjoint_generator);
    SpectralStateOps::normalize_energy(partition_state);
    SpectralStateOps::normalize_energy(partition_tangent_state);
    const LocalCriticalDerivativeLedgerReport critical_derivative_ledger =
        LocalCriticalDerivativeLedger::evaluate(
            active_dynamics, active_objective, partition_state,
            adjoint_viscosity);
    const ShiftedCriticalDensityDiagnostic shifted_density_diagnostic =
        ShiftedCriticalDensityAnalyzer::evaluate(
            active_dynamics, active_objective, partition_state,
            adjoint_viscosity);
    const LocalQuarticIdentityReport local_quartic_identity =
        LocalQuarticIdentityLedger::evaluate(
            active_dynamics, partition_state,
            critical_derivative_ledger);
    const LocalQuarticCommutatorReport local_quartic_commutator =
        LocalQuarticCommutator::evaluate(
            active_dynamics, partition_state,
            critical_derivative_ledger);
    const LocalQuarticProjectedResidualReport local_quartic_projected =
        LocalQuarticProjectedResidual::evaluate(
            active_dynamics, partition_state,
            shifted_density_diagnostic, critical_derivative_ledger);
    const LocalQuarticReducedReport local_quartic_reduced =
        LocalQuarticReducedLedger::evaluate(
            local_quartic_commutator, local_quartic_projected,
            critical_derivative_ledger);
    const LocalQuarticClosureTargetReport local_quartic_closure =
        LocalQuarticClosureTarget::evaluate(
            shifted_density_diagnostic, critical_derivative_ledger,
            local_quartic_reduced);
    const LocalQuarticClosureObjective local_closure_objective(
        active_dynamics);
    const LocalQuarticClosureObjectiveValue local_closure_value =
        local_closure_objective.evaluate(partition_state);
    const Real local_closure_value_error = std::abs(
        local_closure_value.constant_ratio -
        local_quartic_closure.required_constant_ratio) /
        std::max(
            1e-30L,
            std::max(
                std::abs(local_closure_value.constant_ratio),
                std::abs(local_quartic_closure
                             .required_constant_ratio)));
    const LocalQuarticShellReport local_quartic_shells =
        LocalQuarticShellLedger::evaluate(
            active_dynamics, partition_state,
            shifted_density_diagnostic, critical_derivative_ledger);
    const LocalQuarticShellEnvelopeReport local_quartic_envelope =
        LocalQuarticShellEnvelope::analyze(
            partition_state, local_quartic_shells);
    const ShiftedCriticalDensityBudget shifted_density_budget =
        ShiftedCriticalDensityBudgetAnalyzer::evaluate(
            shifted_density_diagnostic, critical_derivative_ledger);
    const SpectralIncrement& partition_tangent =
        partition_tangent_state.velocity;
    active_galerkin.set_compute_threads(1);
    const SpectralIncrement serial_partition_advection =
        active_dynamics.advection_direct_partition(
            partition_state, TriadPartition::nonlocal);
    const SpectralIncrement serial_partition_vjp =
        active_dynamics.advection_vjp_direct_partition(
            partition_state, partition_tangent,
            TriadPartition::nonlocal);
    active_galerkin.set_compute_threads(4);
    const SpectralIncrement parallel_partition_advection =
        active_dynamics.advection_direct_partition(
            partition_state, TriadPartition::nonlocal);
    const SpectralIncrement parallel_partition_vjp =
        active_dynamics.advection_vjp_direct_partition(
            partition_state, partition_tangent,
            TriadPartition::nonlocal);
    active_galerkin.set_compute_threads(1);
    const Real partition_parallel_forward_error = increment_relative_error(
        parallel_partition_advection, serial_partition_advection);
    const Real partition_parallel_vjp_error = increment_relative_error(
        parallel_partition_vjp, serial_partition_vjp);
    const TriadLedgerReport partition_ledger =
        TriadLedger::analyze(partition_state);
    const TriadCommutatorReport partition_commutator =
        TriadCommutator::analyze(partition_state);
    const TriadTailEnvelopeReport partition_tail_envelope =
        TriadTailEnvelope::analyze(partition_state);
    const Real partition_ledger_error = std::abs(
        partition_ledger.signed_local -
        active_objective
            .evaluate(partition_state, TriadPartition::local)
            .signed_vortex_stretching) /
        std::max(1e-30L, std::abs(partition_ledger.signed_local));
    Real tail_envelope_signed = 0.0L;
    bool tail_envelope_bounds_ok = true;
    for (std::size_t role = 0;
         role < separated_low_role_count; ++role) {
        tail_envelope_signed +=
            partition_tail_envelope.signed_stretching_by_low_role[role];
        tail_envelope_bounds_ok = tail_envelope_bounds_ok &&
            partition_tail_envelope.maximum_amplitude_bound_ratio[role] <=
                1.0L + 1e-14L &&
            partition_tail_envelope
                    .maximum_normalized_frequency_ratio[role] <=
                1.0L + 1e-14L;
    }
    const Real tail_envelope_partition_error = std::abs(
        tail_envelope_signed - partition_ledger.signed_nonlocal) /
        std::max(1e-30L, std::abs(partition_ledger.signed_nonlocal));
    const bool partition_parallel_ok =
        partition_parallel_forward_error < 1e-14L &&
        partition_parallel_vjp_error < 1e-14L &&
        partition_ledger_error < 1e-14L &&
        partition_commutator.pairs > 0 &&
        partition_commutator
                .relative_unweighted_cancellation_residual < 1e-14L &&
        partition_commutator.relative_weighted_identity_residual < 1e-14L &&
        partition_commutator.maximum_frequency_gain_ratio <=
            1.0L + 1e-14L &&
        tail_envelope_bounds_ok &&
        tail_envelope_partition_error < 1e-14L;
    const SpectralState partition_plus_state =
        active_dynamics.add_increment(
            partition_state, partition_tangent,
            finite_difference_step);
    const SpectralState partition_minus_state =
        active_dynamics.add_increment(
            partition_state, partition_tangent,
            -finite_difference_step);
    const SpectralIncrement local_closure_gradient =
        local_closure_objective.squared_constant_ratio_gradient(
            partition_state);
    const Real local_closure_directional_adjoint =
        increment_inner_product(
            local_closure_gradient, partition_tangent);
    const Real local_closure_directional_finite_difference =
        (local_closure_objective.evaluate(partition_plus_state)
             .squared_constant_ratio -
         local_closure_objective.evaluate(partition_minus_state)
             .squared_constant_ratio) /
        (2.0L * finite_difference_step);
    const Real local_closure_gradient_error = std::abs(
        local_closure_directional_adjoint -
        local_closure_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(
                std::abs(local_closure_directional_adjoint),
                std::abs(
                    local_closure_directional_finite_difference)));
    const SpectralIncrement signed_closure_gradient =
        local_closure_objective.signed_constant_ratio_gradient(
            partition_state);
    const Real signed_closure_directional_adjoint =
        increment_inner_product(
            signed_closure_gradient, partition_tangent);
    const Real signed_closure_directional_finite_difference =
        (local_closure_objective.evaluate(partition_plus_state)
             .signed_constant_ratio -
         local_closure_objective.evaluate(partition_minus_state)
             .signed_constant_ratio) /
        (2.0L * finite_difference_step);
    const Real signed_closure_gradient_error = std::abs(
        signed_closure_directional_adjoint -
        signed_closure_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(
                std::abs(signed_closure_directional_adjoint),
                std::abs(
                    signed_closure_directional_finite_difference)));
    const SpectralIncrement local_sld_gradient =
        local_closure_objective.signed_local_sld_ratio_gradient(
            partition_state);
    const Real local_sld_directional_adjoint =
        increment_inner_product(local_sld_gradient, partition_tangent);
    const Real local_sld_directional_finite_difference =
        (local_closure_objective.evaluate(partition_plus_state)
             .signed_local_sld_ratio -
         local_closure_objective.evaluate(partition_minus_state)
             .signed_local_sld_ratio) /
        (2.0L * finite_difference_step);
    const Real local_sld_gradient_error = std::abs(
        local_sld_directional_adjoint -
        local_sld_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(
                std::abs(local_sld_directional_adjoint),
                std::abs(local_sld_directional_finite_difference)));
    const TriadSelection doubling_selection =
        TriadSelection::local_equal_low_doubling();
    const LocalSldBlockObjective selected_block_objective(
        active_dynamics, doubling_selection,
        LocalSldBlock::selected_closed);
    const LocalSldBlockObjective complement_block_objective(
        active_dynamics, doubling_selection,
        LocalSldBlock::complement_closed);
    const LocalSldBlockObjective mixed_block_objective(
        active_dynamics, doubling_selection,
        LocalSldBlock::mixed);
    const LocalSldBlockObjective full_block_objective(
        active_dynamics, doubling_selection,
        LocalSldBlock::full);
    auto block_gradient_error = [&](const LocalSldBlockObjective& objective) {
        const SpectralIncrement gradient = objective.gradient(
            partition_state);
        const Real directional_adjoint = increment_inner_product(
            gradient, partition_tangent);
        const Real directional_finite_difference =
            (objective.evaluate(partition_plus_state).block_sld_ratio -
             objective.evaluate(partition_minus_state).block_sld_ratio) /
            (2.0L * finite_difference_step);
        return std::abs(
            directional_adjoint - directional_finite_difference) /
            std::max(
                1e-30L,
                std::max(
                    std::abs(directional_adjoint),
                    std::abs(directional_finite_difference)));
    };
    const Real selected_block_gradient_error = block_gradient_error(
        selected_block_objective);
    const Real complement_block_gradient_error = block_gradient_error(
        complement_block_objective);
    const Real mixed_block_gradient_error = block_gradient_error(
        mixed_block_objective);
    const Real full_block_gradient_error = block_gradient_error(
        full_block_objective);
    const LocalSldBlockObjectiveValue selected_block_value =
        selected_block_objective.evaluate(partition_state);
    const LocalSldBlockObjectiveValue complement_block_value =
        complement_block_objective.evaluate(partition_state);
    const LocalSldBlockObjectiveValue mixed_block_value =
        mixed_block_objective.evaluate(partition_state);
    const Real block_ratio_reconstruction_error = std::abs(
        local_closure_value.signed_local_sld_ratio -
        selected_block_value.block_sld_ratio -
        complement_block_value.block_sld_ratio -
        mixed_block_value.block_sld_ratio) /
        std::max(
            std::abs(local_closure_value.signed_local_sld_ratio),
            1e-30L);
    const LocalSldTrajectoryAdjoint frozen_sld_adjoint(
        active_dynamics, TriadPartition::local);
    const SpectralIncrement frozen_static_gradient =
        local_closure_objective.frozen_signed_local_sld_ratio_gradient(
            partition_state, local_closure_value.initial_frequency,
            local_closure_value.initial_ep_shift);
    const Real frozen_static_directional_adjoint = increment_inner_product(
        frozen_static_gradient, partition_tangent);
    const Real frozen_static_directional_finite_difference =
        (local_closure_objective.frozen_signed_local_sld_ratio(
             partition_plus_state,
             local_closure_value.initial_frequency,
             local_closure_value.initial_ep_shift) -
         local_closure_objective.frozen_signed_local_sld_ratio(
             partition_minus_state,
             local_closure_value.initial_frequency,
             local_closure_value.initial_ep_shift)) /
        (2.0L * finite_difference_step);
    const Real frozen_static_gradient_error = std::abs(
        frozen_static_directional_adjoint -
        frozen_static_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(
                std::abs(frozen_static_directional_adjoint),
                std::abs(frozen_static_directional_finite_difference)));
    constexpr int frozen_sld_steps = 2;
    constexpr Real frozen_sld_dt = 5e-4L;
    const QTrajectoryGradient frozen_sld_gradient =
        frozen_sld_adjoint.terminal_gradient(
            partition_state, adjoint_viscosity,
            frozen_sld_dt, frozen_sld_steps);
    const Real frozen_sld_directional_adjoint = increment_inner_product(
        frozen_sld_gradient.initial_gradient, partition_tangent);
    const Real frozen_sld_directional_finite_difference =
        (frozen_sld_adjoint.terminal_value(
             partition_plus_state, adjoint_viscosity,
             frozen_sld_dt, frozen_sld_steps).terminal_ratio -
         frozen_sld_adjoint.terminal_value(
             partition_minus_state, adjoint_viscosity,
             frozen_sld_dt, frozen_sld_steps).terminal_ratio) /
        (2.0L * finite_difference_step);
    const Real frozen_sld_gradient_error = std::abs(
        frozen_sld_directional_adjoint -
        frozen_sld_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(
                std::abs(frozen_sld_directional_adjoint),
                std::abs(frozen_sld_directional_finite_difference)));
    const QTrajectoryGradient frozen_sld_maximum_gradient =
        frozen_sld_adjoint.maximum_gradient(
            partition_state, adjoint_viscosity,
            frozen_sld_dt, frozen_sld_steps);
    const Real frozen_sld_maximum_directional_adjoint =
        increment_inner_product(
            frozen_sld_maximum_gradient.initial_gradient,
            partition_tangent);
    const Real frozen_sld_maximum_directional_finite_difference =
        (frozen_sld_adjoint.maximum_value(
             partition_plus_state, adjoint_viscosity,
             frozen_sld_dt, frozen_sld_steps).terminal_ratio -
         frozen_sld_adjoint.maximum_value(
             partition_minus_state, adjoint_viscosity,
             frozen_sld_dt, frozen_sld_steps).terminal_ratio) /
        (2.0L * finite_difference_step);
    const Real frozen_sld_maximum_gradient_error = std::abs(
        frozen_sld_maximum_directional_adjoint -
        frozen_sld_maximum_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(
                std::abs(frozen_sld_maximum_directional_adjoint),
                std::abs(
                    frozen_sld_maximum_directional_finite_difference)));
    const QTrajectoryGradient frozen_sld_zero_step =
        frozen_sld_adjoint.terminal_gradient(
            partition_state, adjoint_viscosity,
            frozen_sld_dt, 0);
    const Real frozen_sld_zero_step_gradient_error =
        increment_relative_error(
            frozen_sld_zero_step.initial_gradient,
            local_sld_gradient);
    auto partition_integral_gradient_error = [&](TriadSelection selection) {
        const QTrajectoryGradient partition_gradient =
            active_adjoint.critical_integral_gradient(
                partition_state, adjoint_viscosity, adjoint_dt,
                trajectory_steps, selection);
        const Real directional_adjoint = increment_inner_product(
            partition_gradient.initial_gradient, partition_tangent);
        auto partition_integral = [&](SpectralState state) {
            StaticObjective previous =
                active_objective.evaluate(state, selection);
            Real integral = 0.0L;
            for (int step = 0; step < trajectory_steps; ++step) {
                active_dynamics.rk4_step(
                    state, adjoint_viscosity, adjoint_dt);
                const StaticObjective current =
                    active_objective.evaluate(state, selection);
                integral += 0.5L * adjoint_dt *
                            (previous.critical_integrand +
                             current.critical_integrand);
                previous = current;
            }
            return integral;
        };
        const Real directional_finite_difference =
            (partition_integral(partition_plus_state) -
             partition_integral(partition_minus_state)) /
            (2.0L * finite_difference_step);
        return std::abs(directional_adjoint -
                        directional_finite_difference) /
               std::max(
                   1e-30L,
                   std::max(std::abs(directional_adjoint),
                            std::abs(directional_finite_difference)));
    };
    const Real local_integral_gradient_error =
        partition_integral_gradient_error(TriadPartition::local);
    const Real nonlocal_integral_gradient_error =
        partition_integral_gradient_error(TriadPartition::nonlocal);
    const Real near_nonlocal_integral_gradient_error =
        partition_integral_gradient_error(TriadPartition::near_nonlocal);
    const Real far_nonlocal_integral_gradient_error =
        partition_integral_gradient_error(TriadPartition::far_nonlocal);
    const Real configurable_tail_integral_gradient_error =
        partition_integral_gradient_error(
            TriadSelection::dyadic_tail(1));
    const QTrajectoryGradient local_increase_gradient =
        active_adjoint.critical_increase_gradient(
            partition_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps, TriadPartition::local);
    const Real local_increase_directional = increment_inner_product(
        local_increase_gradient.initial_gradient, partition_tangent);
    auto local_critical_increase = [&](SpectralState state) {
        const Real initial = active_objective
            .evaluate(state, TriadPartition::local)
            .critical_integrand;
        for (int step = 0; step < trajectory_steps; ++step) {
            active_dynamics.rk4_step(
                state, adjoint_viscosity, adjoint_dt);
        }
        return active_objective
                   .evaluate(state, TriadPartition::local)
                   .critical_integrand -
               initial;
    };
    const Real local_increase_finite_difference =
        (local_critical_increase(partition_plus_state) -
         local_critical_increase(partition_minus_state)) /
        (2.0L * finite_difference_step);
    const Real local_increase_gradient_error = std::abs(
        local_increase_directional -
        local_increase_finite_difference) /
        std::max(
            1e-30L,
            std::max(std::abs(local_increase_directional),
                     std::abs(local_increase_finite_difference)));
    constexpr Real local_log_gain_shift = 1e-6L;
    const QTrajectoryGradient local_log_gain_gradient =
        active_adjoint.critical_log_gain_gradient(
            partition_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps, TriadPartition::local,
            local_log_gain_shift);
    const Real local_log_gain_directional = increment_inner_product(
        local_log_gain_gradient.initial_gradient, partition_tangent);
    auto local_critical_log_gain = [&](SpectralState state) {
        const Real initial = active_objective
            .evaluate(state, TriadPartition::local)
            .critical_integrand;
        for (int step = 0; step < trajectory_steps; ++step) {
            active_dynamics.rk4_step(
                state, adjoint_viscosity, adjoint_dt);
        }
        const Real terminal = active_objective
            .evaluate(state, TriadPartition::local)
            .critical_integrand;
        return std::log(
            (terminal + local_log_gain_shift) /
            (initial + local_log_gain_shift));
    };
    const Real local_log_gain_finite_difference =
        (local_critical_log_gain(partition_plus_state) -
         local_critical_log_gain(partition_minus_state)) /
        (2.0L * finite_difference_step);
    const Real local_log_gain_gradient_error = std::abs(
        local_log_gain_directional -
        local_log_gain_finite_difference) /
        std::max(
            1e-30L,
            std::max(std::abs(local_log_gain_directional),
                     std::abs(local_log_gain_finite_difference)));
    const QTrajectoryGradient local_ep_log_gain_gradient =
        active_adjoint.critical_ep_log_gain_gradient(
            partition_state, adjoint_viscosity, adjoint_dt,
            trajectory_steps, TriadPartition::local);
    const Real local_ep_log_gain_directional = increment_inner_product(
        local_ep_log_gain_gradient.initial_gradient, partition_tangent);
    auto local_critical_ep_log_gain = [&](SpectralState state) {
        const StaticObjective initial = active_objective.evaluate(
            state, TriadPartition::local);
        const Real shift = initial.energy * initial.palinstrophy;
        for (int step = 0; step < trajectory_steps; ++step) {
            active_dynamics.rk4_step(
                state, adjoint_viscosity, adjoint_dt);
        }
        const Real terminal = active_objective
            .evaluate(state, TriadPartition::local)
            .critical_integrand;
        return std::log1p(
            (terminal - initial.critical_integrand) /
            (initial.critical_integrand + shift));
    };
    const Real local_ep_log_gain_finite_difference =
        (local_critical_ep_log_gain(partition_plus_state) -
         local_critical_ep_log_gain(partition_minus_state)) /
        (2.0L * finite_difference_step);
    const Real local_ep_log_gain_gradient_error = std::abs(
        local_ep_log_gain_directional -
        local_ep_log_gain_finite_difference) /
        std::max(
            1e-30L,
            std::max(std::abs(local_ep_log_gain_directional),
                     std::abs(local_ep_log_gain_finite_difference)));
    const EvolutionResult configurable_tail_evolution =
        active_trajectory_analyzer.evolve(
            partition_state, adjoint_viscosity,
            adjoint_dt * static_cast<Real>(trajectory_steps),
            adjoint_dt, true, 1);
    const Real configurable_tail_trajectory_error = std::abs(
        configurable_tail_evolution.integral_selected_gap_tail_critical -
        configurable_tail_evolution.integral_nonlocal_critical) /
        std::max(
            1e-30L,
            std::abs(configurable_tail_evolution
                         .integral_nonlocal_critical));
    SpectralState far_partition_state =
        SpectralStateFactory::random(3, adjoint_generator);
    SpectralState far_partition_tangent_state =
        SpectralStateFactory::random(3, adjoint_generator);
    SpectralStateOps::normalize_energy(far_partition_state);
    SpectralStateOps::normalize_energy(far_partition_tangent_state);
    const SpectralIncrement far_partition_gradient =
        active_objective.critical_integrand_gradient(
            far_partition_state, TriadPartition::far_nonlocal);
    const Real far_partition_directional_gradient = increment_inner_product(
        far_partition_gradient, far_partition_tangent_state.velocity);
    const SpectralState far_partition_plus = active_dynamics.add_increment(
        far_partition_state, far_partition_tangent_state.velocity,
        finite_difference_step);
    const SpectralState far_partition_minus = active_dynamics.add_increment(
        far_partition_state, far_partition_tangent_state.velocity,
        -finite_difference_step);
    const Real far_partition_directional_finite_difference =
        (active_objective
             .evaluate(far_partition_plus, TriadPartition::far_nonlocal)
             .critical_integrand -
         active_objective
             .evaluate(far_partition_minus, TriadPartition::far_nonlocal)
             .critical_integrand) /
        (2.0L * finite_difference_step);
    const Real far_partition_static_gradient_error = std::abs(
        far_partition_directional_gradient -
        far_partition_directional_finite_difference) /
        std::max(
            1e-30L,
            std::max(std::abs(far_partition_directional_gradient),
                     std::abs(
                         far_partition_directional_finite_difference)));
    const bool triad_partition_ok =
        TriadPartitioner::is_local(
            WaveVector{1, 0, 0}, WaveVector{0, 1, 0},
            WaveVector{1, 1, 0}) &&
        !TriadPartitioner::is_local(
            WaveVector{1, 0, 0}, WaveVector{2, 0, 0},
            WaveVector{3, 0, 0}) &&
        TriadPartitioner::dyadic_gap(
            WaveVector{1, 0, 0}, WaveVector{2, 0, 0},
            WaveVector{3, 0, 0}) == 1 &&
        TriadPartitioner::dyadic_gap(
            WaveVector{1, 0, 0}, WaveVector{7, 0, 0},
            WaveVector{8, 0, 0}) == 2 &&
        TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{7, 0, 0},
            WaveVector{8, 0, 0}, TriadPartition::far_nonlocal) &&
        TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{2, 0, 0},
            WaveVector{3, 0, 0}, TriadPartition::nonlocal) &&
        TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{9, 0, 0},
            WaveVector{10, 0, 0}, TriadSelection::dyadic_tail(3)) &&
        !TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{7, 0, 0},
            WaveVector{8, 0, 0}, TriadSelection::dyadic_tail(3));
    const bool partition_integral_gradients_ok =
        triad_partition_ok && partition_parallel_ok &&
        critical_derivative_ledger.finite &&
        local_quartic_identity.finite && local_quartic_shells.finite &&
        local_quartic_commutator.finite &&
        local_quartic_projected.finite &&
        local_quartic_reduced.finite &&
        local_quartic_closure.finite &&
        local_closure_value.finite &&
        local_closure_value.factorization_relative_error < 1e-12L &&
        local_closure_value_error < 1e-12L &&
        local_closure_gradient_error < 1e-9L &&
        signed_closure_gradient_error < 1e-9L &&
        local_sld_gradient_error < 1e-9L &&
        selected_block_gradient_error < 1e-9L &&
        complement_block_gradient_error < 1e-9L &&
        mixed_block_gradient_error < 1e-9L &&
        full_block_gradient_error < 1e-9L &&
        block_ratio_reconstruction_error < 1e-12L &&
        frozen_sld_gradient_error < 1e-9L &&
        frozen_sld_maximum_gradient_error < 1e-9L &&
        frozen_static_gradient_error < 1e-9L &&
        frozen_sld_zero_step_gradient_error < 1e-12L &&
        local_quartic_envelope.all_bounds_hold &&
        shifted_density_budget.finite &&
        local_integral_gradient_error < 1e-9L &&
        nonlocal_integral_gradient_error < 1e-9L &&
        near_nonlocal_integral_gradient_error < 1e-9L &&
        far_nonlocal_integral_gradient_error < 1e-9L &&
        configurable_tail_integral_gradient_error < 1e-9L &&
        local_increase_gradient_error < 1e-9L &&
        local_log_gain_gradient_error < 1e-9L &&
        local_ep_log_gain_gradient_error < 1e-9L &&
        configurable_tail_trajectory_error < 1e-14L &&
        far_partition_static_gradient_error < 1e-9L;
    GradientSearchOptions gradient_options;
    gradient_options.iterations = 3;
    gradient_options.line_search_steps = 8;
    gradient_options.trajectory_steps = trajectory_steps;
    gradient_options.viscosity = adjoint_viscosity;
    gradient_options.time_step = adjoint_dt;
    gradient_options.initial_step = 0.2L;
    gradient_options.method = "lbfgs";
    const GradientSearchResult gradient_search =
        active_gradient_adversary.maximize_q(
            adjoint_state, gradient_options);
    const Real gradient_energy_error = std::abs(
        SpectralStateOps::energy(gradient_search.state) -
        SpectralStateOps::energy(adjoint_state));
    Real gradient_constraint_error = 0.0L;
    for (std::size_t mode = 0;
         mode < gradient_search.state.waves.size(); ++mode) {
        const WaveVector wave = gradient_search.state.waves[mode];
        gradient_constraint_error = std::max(
            gradient_constraint_error,
            std::abs(wave_dot(
                wave, gradient_search.state.velocity[mode])));
        const std::size_t negative =
            gradient_search.state.index.at(-wave);
        for (std::size_t component = 0; component < 3; ++component) {
            gradient_constraint_error = std::max(
                gradient_constraint_error,
                std::abs(gradient_search.state.velocity[negative]
                             [component] -
                         std::conj(gradient_search.state.velocity[mode]
                                      [component])));
        }
    }
    const bool gradient_search_ok =
        gradient_search.objective >= gradient_search.initial_objective &&
        gradient_search.accepted_steps > 0 &&
        gradient_energy_error < 1e-14L &&
        gradient_constraint_error < 1e-15L;
    GradientSearchOptions closure_gradient_options;
    closure_gradient_options.iterations = 3;
    closure_gradient_options.line_search_steps = 12;
    closure_gradient_options.trajectory_steps = 0;
    closure_gradient_options.initial_step = 0.05L;
    closure_gradient_options.objective = "local-closure-ratio";
    closure_gradient_options.method = "lbfgs";
    const GradientSearchResult closure_gradient_search =
        active_gradient_adversary.maximize_q(
            partition_state, closure_gradient_options);
    const bool closure_gradient_search_ok =
        closure_gradient_search.objective >=
            closure_gradient_search.initial_objective &&
        closure_gradient_search.accepted_steps > 0 &&
        std::isfinite(closure_gradient_search.objective);
    closure_gradient_options.objective =
        "local-signed-closure-ratio";
    const GradientSearchResult signed_closure_gradient_search =
        active_gradient_adversary.maximize_q(
            partition_state, closure_gradient_options);
    const bool signed_closure_gradient_search_ok =
        signed_closure_gradient_search.objective >=
            signed_closure_gradient_search.initial_objective &&
        signed_closure_gradient_search.accepted_steps > 0 &&
        std::isfinite(signed_closure_gradient_search.objective);
    closure_gradient_options.objective = "local-sld-ratio";
    const GradientSearchResult local_sld_gradient_search =
        active_gradient_adversary.maximize_q(
            partition_state, closure_gradient_options);
    const bool local_sld_gradient_search_ok =
        local_sld_gradient_search.objective >=
            local_sld_gradient_search.initial_objective &&
        local_sld_gradient_search.accepted_steps > 0 &&
        std::isfinite(local_sld_gradient_search.objective);
    LocalSldCyclicAnsatzOptions cyclic_options;
    cyclic_options.coarse_samples = 512;
    cyclic_options.refinement_iterations = 48;
    const LocalSldCyclicAnsatzReport cyclic_ansatz =
        LocalSldCyclicAnsatz::optimize(cyclic_options);
    const bool cyclic_ansatz_ok = cyclic_ansatz.value.finite &&
        cyclic_ansatz.value.signed_local_sld_ratio > 7.7e-4L &&
        std::abs(cyclic_ansatz.basis_inner_product) < 1e-14L &&
        cyclic_ansatz.pure_axis_identity_error < 1e-14L;
    LocalSldCyclicTrajectoryOptions cyclic_trajectory_options;
    cyclic_trajectory_options.coarse_samples = 32;
    cyclic_trajectory_options.refinement_iterations = 2;
    cyclic_trajectory_options.trajectory_steps = 2;
    cyclic_trajectory_options.threads = 2;
    const LocalSldCyclicTrajectoryReport cyclic_trajectory_ansatz =
        LocalSldCyclicTrajectoryAnsatz::optimize(
            cyclic_trajectory_options);
    const bool cyclic_trajectory_ansatz_ok =
        cyclic_trajectory_ansatz.value.finite &&
        cyclic_trajectory_ansatz.refined_value.finite &&
        cyclic_trajectory_ansatz.value.steps >= 0 &&
        cyclic_trajectory_ansatz.value.steps <= 2 &&
        std::isfinite(cyclic_trajectory_ansatz.restricted_gradient) &&
        std::isfinite(
            cyclic_trajectory_ansatz.projected_full_gradient_norm) &&
        cyclic_trajectory_ansatz.time_step_relative_error < 1e-4L;
    LocalSldCyclicKrylovOptions cyclic_krylov_options;
    cyclic_krylov_options.cutoff = 2;
    cyclic_krylov_options.warm_angle_samples = 16;
    cyclic_krylov_options.iterations = 1;
    cyclic_krylov_options.line_search_steps = 2;
    cyclic_krylov_options.trajectory_steps = 2;
    cyclic_krylov_options.threads = 2;
    cyclic_krylov_options.backend = "direct";
    const LocalSldCyclicKrylovReport cyclic_krylov_ansatz =
        LocalSldCyclicKrylovAnsatz::optimize(cyclic_krylov_options);
    const bool cyclic_krylov_ansatz_ok =
        cyclic_krylov_ansatz.value.finite &&
        cyclic_krylov_ansatz.refined_value.finite &&
        cyclic_krylov_ansatz.maximum_gram_error < 1e-14L &&
        std::isfinite(cyclic_krylov_ansatz.restricted_gradient_norm) &&
        std::isfinite(
            cyclic_krylov_ansatz.projected_full_gradient_norm) &&
        cyclic_krylov_ansatz.time_step_relative_error < 1e-4L;
    const LocalSldResponseHierarchyReport response_hierarchy =
        LocalSldResponseHierarchy::analyze(
            active_dynamics, cyclic_ansatz.state, 4);
    const bool response_hierarchy_ok =
        response_hierarchy.constructed_depth == 4 &&
        response_hierarchy.maximum_gram_error < 1e-14L &&
        response_hierarchy.final_projection_energy > 0.999999L &&
        response_hierarchy.final_projection_energy < 1.000001L;
    const LocalSldResponseFamilyReport response_family =
        LocalSldResponseFamily::analyze(
            active_dynamics,
            {{"reference-a", cyclic_ansatz.state},
             {"reference-b", cyclic_ansatz.state}},
            4, false, false);
    const bool response_family_ok =
        response_family.rows.size() == 2 &&
        response_family.depth == 4 &&
        response_family.maximum_coefficient_l2_difference < 1e-18L &&
        std::abs(
            response_family.rows[0].final_projection_energy -
            response_family.rows[1].final_projection_energy) < 1e-18L &&
        response_family.finite_cutoff_family_is_not_a_proof;
    const LocalSldTrajectoryEvaluatorReport trajectory_evaluation =
        LocalSldTrajectoryEvaluator::evaluate(
            active_dynamics, cyclic_ansatz.state,
            0.1L, 0.001L, 2);
    const bool trajectory_evaluation_ok =
        trajectory_evaluation.finite &&
        trajectory_evaluation.maximum.terminal_ratio >=
            trajectory_evaluation.initial.terminal_ratio &&
        trajectory_evaluation.maximum.terminal_ratio >=
            trajectory_evaluation.terminal.terminal_ratio &&
        trajectory_evaluation.time_step_relative_error < 1e-4L;
    const LocalSldSignatureBlockReport signature_block =
        LocalSldSignatureBlock::analyze(
            active_dynamics, cyclic_ansatz.state, {1, 1, 2});
    const bool signature_block_ok = signature_block.finite &&
        signature_block.dominant_interactions > 0 &&
        signature_block.quotient_reconstruction_error < 1e-14L &&
        TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{0, 1, 0},
            WaveVector{1, 1, 0},
            TriadSelection::local_signature(1, 1, 2)) &&
        !TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{0, 1, 0},
            WaveVector{1, 1, 0},
            TriadSelection::local_without_signature(1, 1, 2)) &&
        TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{0, 1, 0},
            WaveVector{1, 1, 0},
            TriadSelection::local_equal_low_doubling()) &&
        !TriadPartitioner::includes(
            WaveVector{1, 0, 0}, WaveVector{0, 1, 0},
            WaveVector{1, 1, 0},
            TriadSelection::local_without_equal_low_doubling());
    const AdversaryResult adversary =
        optimize_static_depletion(1, 1, 2, 0.1L, 11);
    const bool adversary_ok = adversary.modes == 26 &&
                              std::abs(adversary.objective.energy - 1.0L) < 1e-15L &&
                              std::isfinite(adversary.objective.energy_level_quantity) &&
                              adversary.objective.energy_level_quantity >= 0.0L;
    const EvolutionResult evolution =
        active_trajectory_analyzer.evolve(adversary.state, 0.1L, 0.002L, 0.001L);
    DynamicAdversaryOptions dynamic_class_options;
    dynamic_class_options.generations = 0;
    dynamic_class_options.viscosity = 0.1L;
    dynamic_class_options.final_time = 0.002L;
    dynamic_class_options.time_step = 0.001L;
    dynamic_class_options.objective = "critical-integral";
    dynamic_class_options.seed = 17;
    const DynamicAdversaryEnsemble dynamic_class_adversary("direct", 2);
    const DynamicAdversaryResult dynamic_class_result =
        dynamic_class_adversary.optimize(
            adversary.state, nullptr, dynamic_class_options, 2);
    const bool dynamic_class_ok =
        dynamic_class_result.refined_evolution.finite &&
        dynamic_class_result.restart_objectives.size() == 2 &&
        dynamic_class_result.winning_restart >= 0 &&
        dynamic_class_result.winning_restart < 2 &&
        std::isfinite(dynamic_class_result.search_final_objective) &&
        dynamic_class_result.time_step_relative_error < 1e-4L;
    const QDerivativeDiagnostic q_derivative =
        active_trajectory_analyzer.evaluate_q_derivative(adversary.state, 0.1L);
    const bool q_derivative_ok = q_derivative.valid &&
        q_derivative.relative_refinement_error < 1e-6L;
    const bool evolution_ok = evolution.finite && evolution.steps == 2 &&
                              evolution.final_energy <= evolution.initial_energy &&
                              std::abs(evolution.energy_balance_residual) < 1e-10L;
    out << "rational/scaling test: " << (rational_ok && scaling_ok ? "PASS" : "FAIL")
        << " (minimum gamma=" << scaling.minimum_young_power.str() << ")\n"
        << "concentration scaling test: "
        << (concentration_ok ? "PASS" : "FAIL")
        << " (Q exponent=" << concentration.fixed_energy_pointwise_q.str()
        << ", integral exponent=" << concentration.natural_integrated_l4.str()
        << ")\n"
        << "strong L4 reduction test: " << (strong_l4_ok ? "PASS" : "FAIL")
        << " (integral D4Z2 <= sup(Q)*E0/(2nu))\n"
        << "shifted local density scaling test: "
        << (shifted_density_ok ? "PASS" : "FAIL")
        << " (C and E0P0 amplitude="
        << shifted_density.density_amplitude_degree.str()
        << ", scaling=" << shifted_density.density_scaling_exponent.str()
        << ", conditional energy closure)\n"
        << "dyadic tail scaling test: "
        << (dyadic_tail_scaling_ok ? "PASS" : "FAIL")
        << " (gap decay=" << dyadic_tail.low_advecting_gap_decay.str()
        << ", remaining Z power="
        << dyadic_tail.l4_density_enstrophy_power.str()
        << ", energy closure="
        << (dyadic_tail.energy_identity_closes_time_integral
                ? "YES"
                : "no")
        << ", moving-gap remainder Z^"
        << dyadic_tail.moving_gap_remaining_enstrophy_power.str()
        << ")\n"
        << "dyadic shell sequence test: "
        << (dyadic_shell_bounds_ok ? "PASS" : "FAIL")
        << " (high="
        << static_cast<double>(dyadic_shell_bounds.maximum_high_moment_ratio)
        << ", one-gain="
        << static_cast<double>(dyadic_shell_bounds.maximum_one_gain_tail_ratio)
        << ", three-gain="
        << static_cast<double>(dyadic_shell_bounds.maximum_three_gain_tail_ratio)
        << ")\n"
        << "periodic shell geometry test: "
        << (periodic_shell_geometry_ok ? "PASS" : "FAIL")
        << " (count="
        << static_cast<double>(periodic_shell_geometry.maximum_count_ratio)
        << ", overlap2="
        << static_cast<double>(
               periodic_shell_geometry.maximum_one_gain_overlap_ratio)
        << ", overlap1="
        << static_cast<double>(
               periodic_shell_geometry.maximum_three_gain_overlap_ratio)
        << ", C1="
        << static_cast<double>(periodic_shell_geometry.ft1_one_gain_constant)
        << ")\n"
        << "explicit periodic FT-1 test: "
        << (periodic_tail_bound_ok ? "PASS" : "FAIL")
        << " (cutoff=" << periodic_tail_bound.cutoff
        << ", samples=" << periodic_tail_bound.samples
        << ", max ratio="
        << static_cast<double>(periodic_tail_bound.maximum_bound_ratio)
        << ")\n"
        << "moving far-tail closure test: "
        << (far_tail_closure_ok ? "PASS" : "FAIL")
        << " (samples=" << far_tail_closure.samples
        << ", max remainder ratio="
        << static_cast<double>(
               far_tail_closure.maximum_normalized_remainder_ratio)
        << ")\n"
        << "transition block scaling test: "
        << (transition_block_scaling_ok ? "PASS" : "FAIL")
        << " (remainder=log(1+Z)^"
        << transition_block_scaling.post_young_logarithm_power.str()
        << " Z^"
        << transition_block_scaling.post_young_enstrophy_power.str()
        << ", required depletion=Z^(-"
        << transition_block_scaling.required_pointwise_depletion_power.str()
        << "))\n"
        << "moving gap controller test: "
        << (moving_gap_controller_ok ? "PASS" : "FAIL")
        << " (m(1.1)="
        << MovingGapController::decide(1.1L, 2).minimum_gap
        << ", m(4.1)="
        << MovingGapController::decide(4.1L, 2).minimum_gap
        << ")\n"
        << "spectral triad test: " << (triad_ok ? "PASS" : "FAIL")
        << " (energy residual="
        << static_cast<double>(triads.maximum_normalized_energy_residual) << ")\n"
        << "helical triad ledger test: "
        << (helical_ok ? "PASS" : "FAIL")
        << " (velocity="
        << static_cast<double>(
               helical.relative_velocity_reconstruction_residual)
        << ", total="
        << static_cast<double>(helical.relative_total_reconstruction_residual)
        << ", local="
        << static_cast<double>(helical.relative_local_reconstruction_residual)
        << ")\n"
        << "helical gap ledger test: "
        << (helical_gap_ok ? "PASS" : "FAIL")
        << " (gap="
        << static_cast<double>(
               helical_gaps.maximum_gap_reconstruction_residual)
        << ", total="
        << static_cast<double>(helical_gaps.total_reconstruction_residual)
        << ")\n"
        << "local triad symmetrizer test: "
        << (local_symmetry_ok ? "PASS" : "FAIL")
        << " (energy="
        << static_cast<double>(
               local_symmetry.maximum_energy_cancellation_residual)
        << ", reconstruction="
        << static_cast<double>(local_symmetry.local_reconstruction_residual)
        << ", spread="
        << static_cast<double>(
               local_symmetry.maximum_frequency_spread_bound_ratio)
        << ", effective signatures="
        << static_cast<double>(
               local_symmetry.effective_coherent_signature_count)
        << ", signed amplification="
        << static_cast<double>(
               local_symmetry.signed_signature_amplification)
        << ")\n"
        << "orthogonal triad closure test: "
        << (orthogonal_geometry_ok ? "PASS" : "FAIL")
        << " (input degree="
        << static_cast<double>(
               orthogonal_geometry.maximum_input_degree_ratio)
        << ", target degree="
        << static_cast<double>(
               orthogonal_geometry.maximum_target_degree_ratio)
        << ", high-frequency power="
        << orthogonal_closure.transfer_to_viscosity_frequency_power.str()
        << ")\n"
        << "local signature closure test: "
        << (local_signature_geometry_ok ? "PASS" : "FAIL")
        << " (input degree="
        << static_cast<double>(
               local_signature_geometry.maximum_input_degree_ratio)
        << ", target degree="
        << static_cast<double>(
               local_signature_geometry.maximum_target_degree_ratio)
        << ", finite="
        << signature_family_closure.finite_signature_family
               .transfer_frequency_power.str()
        << ", critical="
        << signature_family_closure.critical_signature_family
               .transfer_frequency_power.str()
        << ", dense="
        << signature_family_closure.dense_signature_family
               .transfer_frequency_power.str()
        << ")\n"
        << "local signature objective gradient test: "
        << (local_signature_objective_ok ? "PASS" : "FAIL")
        << " (ledger="
        << static_cast<double>(local_signature_objective_error)
        << ", gradient="
        << static_cast<double>(local_signature_gradient_error)
        << ", transfer gradient="
        << static_cast<double>(local_signature_transfer_gradient_error)
        << ", factorization="
        << static_cast<double>(
               local_signature_density.factorization_residual)
        << ", A_sig="
        << static_cast<double>(
               local_signature_objective.signed_amplification)
        << ")\n"
        << "pure helical local test: "
        << (pure_helical_ok ? "PASS" : "FAIL")
        << " (plus local="
        << static_cast<double>(
               positive_helical.homochiral_local_stretching)
        << ", minus local="
        << static_cast<double>(
               negative_helical.homochiral_local_stretching)
        << ")\n"
        << "helical sector objective gradient test: "
        << (helical_sector_objective_ok ? "PASS" : "FAIL")
        << " (partition="
        << static_cast<double>(helical_sector_partition_error)
        << ", signed="
        << static_cast<double>(helical_signed_gradient_error)
        << ", critical="
        << static_cast<double>(helical_critical_gradient_error)
        << ", broad="
        << static_cast<double>(broad_helical_signed_gradient_error)
        << "/"
        << static_cast<double>(broad_helical_critical_gradient_error)
        << ")\n"
        << "helical sector adversary test: "
        << (helical_adversary_ok ? "PASS" : "FAIL")
        << " (objective="
        << static_cast<double>(helical_adversary.initial_objective)
        << " -> " << static_cast<double>(helical_adversary.objective)
        << ", accepted=" << helical_adversary.accepted_steps << ")\n"
        << "helical trajectory adjoint test: "
        << (helical_trajectory_adjoint_ok ? "PASS" : "FAIL")
        << " (relative error="
        << static_cast<double>(helical_trajectory_gradient_error)
        << ", checkpoints="
        << helical_trajectory_gradient.checkpoint_count << ")\n"
        << "helical trajectory adversary test: "
        << (helical_trajectory_adversary_ok ? "PASS" : "FAIL")
        << " (objective="
        << static_cast<double>(
               helical_trajectory_adversary.initial_objective)
        << " -> "
        << static_cast<double>(helical_trajectory_adversary.objective)
        << ", accepted="
        << helical_trajectory_adversary.accepted_steps << ")\n"
        << "dealiased FFT/direct test: " << (fft_ok ? "PASS" : "FAIL")
        << " (relative error=" << static_cast<double>(fft_relative_error) << ")\n"
        << "FFT adjoint/direct oracle test: "
        << (fft_adjoint_ok ? "PASS" : "FAIL")
        << " (jvp=" << static_cast<double>(fft_jvp_oracle_error)
        << ", vjp=" << static_cast<double>(fft_vjp_oracle_error)
        << ", duality=" << static_cast<double>(fft_adjoint_duality_error)
        << ")\n"
        << "discrete adjoint test: " << (adjoint_ok ? "PASS" : "FAIL")
        << " (rhs jvp=" << static_cast<double>(rhs_jvp_error)
        << ", rhs duality=" << static_cast<double>(rhs_duality_error)
        << ", rk4 jvp=" << static_cast<double>(rk4_jvp_error)
        << ", rk4 duality=" << static_cast<double>(rk4_duality_error)
        << ")\n"
        << "Q objective gradient test: "
        << (q_gradient_ok ? "PASS" : "FAIL")
        << " (relative error=" << static_cast<double>(q_gradient_error)
        << ")\n"
        << "checkpointed trajectory gradient test: "
        << (trajectory_gradient_ok ? "PASS" : "FAIL")
        << " (relative error="
        << static_cast<double>(trajectory_gradient_error)
        << ", checkpoints=" << trajectory_gradient.checkpoint_count << ")\n"
        << "Q-gain trajectory gradient test: "
        << (q_gain_gradient_ok ? "PASS" : "FAIL")
        << " (relative error=" << static_cast<double>(q_gain_gradient_error)
        << ")\n"
        << "Q-increase trajectory gradient test: "
        << (q_increase_gradient_ok ? "PASS" : "FAIL")
        << " (relative error="
        << static_cast<double>(q_increase_gradient_error)
        << ", divergence="
        << static_cast<double>(q_increase_divergence_residual)
        << ", reality="
        << static_cast<double>(q_increase_reality_residual) << ")\n"
        << "critical L4 integral gradient test: "
        << (critical_integral_gradient_ok ? "PASS" : "FAIL")
        << " (relative error="
        << static_cast<double>(critical_integral_gradient_error) << ")\n"
        << "partitioned L4 integral gradient test: "
        << (partition_integral_gradients_ok ? "PASS" : "FAIL")
        << " (classification=" << (triad_partition_ok ? "exact" : "FAILED")
        << ", parallel_forward="
        << static_cast<double>(partition_parallel_forward_error)
        << ", parallel_vjp="
        << static_cast<double>(partition_parallel_vjp_error)
        << ", ledger=" << static_cast<double>(partition_ledger_error)
        << ", commutator="
        << static_cast<double>(partition_commutator
                                   .relative_weighted_identity_residual)
        << ", frequency_gain="
        << static_cast<double>(
               partition_commutator.maximum_frequency_gain_ratio)
        << ", tail_envelope="
        << static_cast<double>(tail_envelope_partition_error)
        << ", local=" << static_cast<double>(local_integral_gradient_error)
        << ", nonlocal="
        << static_cast<double>(nonlocal_integral_gradient_error)
        << ", near="
        << static_cast<double>(near_nonlocal_integral_gradient_error)
        << ", far="
        << static_cast<double>(far_nonlocal_integral_gradient_error)
        << ", configurable_tail="
        << static_cast<double>(configurable_tail_integral_gradient_error)
        << ", local_increase="
        << static_cast<double>(local_increase_gradient_error)
        << ", local_log_gain="
        << static_cast<double>(local_log_gain_gradient_error)
        << ", local_E0P0_log_gain="
        << static_cast<double>(local_ep_log_gain_gradient_error)
        << ", tail_trajectory="
        << static_cast<double>(configurable_tail_trajectory_error)
        << ", far_static="
        << static_cast<double>(far_partition_static_gradient_error)
        << ")\n"
        << "local critical derivative ledger test: "
        << (critical_derivative_ledger.finite ? "PASS" : "FAIL")
        << " (relative reconstruction error="
        << static_cast<double>(critical_derivative_ledger
                                   .relative_reconstruction_error)
        << ", nonlinear="
        << static_cast<double>(critical_derivative_ledger
                 .reconstructed_density_derivative.nonlinear)
        << ", viscous="
        << static_cast<double>(critical_derivative_ledger
                 .reconstructed_density_derivative.viscous)
        << ", normalized budget="
        << static_cast<double>(
               shifted_density_budget.reconstructed_normalized_rate)
        << ", quartic square error="
        << static_cast<double>(local_quartic_identity
                                   .local_outer_negative_square_error)
        << ", shell error="
        << static_cast<double>(local_quartic_shells
                                   .normalized_reconstruction_error)
        << ", envelope ratio="
        << static_cast<double>(
               local_quartic_envelope.maximum_bound_ratio)
        << ", global ZP ratio="
        << static_cast<double>(local_quartic_envelope.global_zp_ratio)
        << ", commutator error="
        << static_cast<double>(local_quartic_commutator
                                   .identity_relative_error)
        << ", projected residual error="
        << static_cast<double>(local_quartic_projected
                                   .normalized_reconstruction_error)
        << ", reduced error="
        << static_cast<double>(local_quartic_reduced
                                   .normalized_reconstruction_error)
        << ", polynomial reduced error="
        << static_cast<double>(local_quartic_reduced
                                   .polynomial_reconstruction_error)
        << ", closure C="
        << static_cast<double>(local_quartic_closure
                                   .required_constant_ratio)
        << ", closure value error="
        << static_cast<double>(local_closure_value_error)
        << ", closure gradient error="
        << static_cast<double>(local_closure_gradient_error)
        << ", signed closure gradient error="
        << static_cast<double>(signed_closure_gradient_error)
        << ", direct SLD gradient error="
        << static_cast<double>(local_sld_gradient_error)
        << ", block gradients="
        << static_cast<double>(selected_block_gradient_error) << '/'
        << static_cast<double>(complement_block_gradient_error) << '/'
        << static_cast<double>(mixed_block_gradient_error) << '/'
        << static_cast<double>(full_block_gradient_error)
        << ", block reconstruction="
        << static_cast<double>(block_ratio_reconstruction_error)
        << ", frozen trajectory gradient="
        << static_cast<double>(frozen_sld_gradient_error)
        << ", frozen maximum="
        << static_cast<double>(frozen_sld_maximum_gradient_error)
        << ", frozen static="
        << static_cast<double>(frozen_static_gradient_error)
        << ", frozen zero-step="
        << static_cast<double>(frozen_sld_zero_step_gradient_error)
        << ", SLD factorization error="
        << static_cast<double>(
               local_closure_value.factorization_relative_error)
        << ")\n"
        << "projected gradient adversary test: "
        << (gradient_search_ok ? "PASS" : "FAIL")
        << " (Q " << static_cast<double>(gradient_search.initial_objective)
        << " -> " << static_cast<double>(gradient_search.objective)
        << ", accepted=" << gradient_search.accepted_steps
        << ", energy error=" << static_cast<double>(gradient_energy_error)
        << ", constraint error="
        << static_cast<double>(gradient_constraint_error)
        << ")\n"
        << "local closure gradient adversary test: "
        << (closure_gradient_search_ok ? "PASS" : "FAIL")
        << " (C "
        << std::sqrt(static_cast<double>(
               closure_gradient_search.initial_objective))
        << " -> "
        << std::sqrt(static_cast<double>(closure_gradient_search.objective))
        << ", accepted=" << closure_gradient_search.accepted_steps
        << ")\n"
        << "signed closure gradient adversary test: "
        << (signed_closure_gradient_search_ok ? "PASS" : "FAIL")
        << " (c "
        << static_cast<double>(
               signed_closure_gradient_search.initial_objective)
        << " -> "
        << static_cast<double>(signed_closure_gradient_search.objective)
        << ", accepted="
        << signed_closure_gradient_search.accepted_steps
        << ")\n"
        << "direct local SLD gradient adversary test: "
        << (local_sld_gradient_search_ok ? "PASS" : "FAIL")
        << " (ratio "
        << static_cast<double>(
               local_sld_gradient_search.initial_objective)
        << " -> "
        << static_cast<double>(local_sld_gradient_search.objective)
        << ", accepted=" << local_sld_gradient_search.accepted_steps
        << ")\n"
        << "local SLD cyclic ansatz test: "
        << (cyclic_ansatz_ok ? "PASS" : "FAIL")
        << " (ratio="
        << static_cast<double>(
               cyclic_ansatz.value.signed_local_sld_ratio)
        << ", basis inner product="
        << static_cast<double>(cyclic_ansatz.basis_inner_product)
        << ", pure axis c="
        << static_cast<double>(
               cyclic_ansatz.pure_axis_value.signed_constant_ratio)
        << ", pure axis error="
        << static_cast<double>(cyclic_ansatz.pure_axis_identity_error)
        << ")\n"
        << "local SLD cyclic trajectory ansatz test: "
        << (cyclic_trajectory_ansatz_ok ? "PASS" : "FAIL")
        << " (ratio="
        << static_cast<double>(
               cyclic_trajectory_ansatz.value.terminal_ratio)
        << ", peak step=" << cyclic_trajectory_ansatz.value.steps
        << ", time refinement="
        << static_cast<double>(
               cyclic_trajectory_ansatz.time_step_relative_error)
        << ")\n"
        << "local SLD cyclic Krylov ansatz test: "
        << (cyclic_krylov_ansatz_ok ? "PASS" : "FAIL")
        << " (ratio="
        << static_cast<double>(
               cyclic_krylov_ansatz.value.terminal_ratio)
        << ", Gram error="
        << static_cast<double>(
               cyclic_krylov_ansatz.maximum_gram_error)
        << ", accepted=" << cyclic_krylov_ansatz.accepted_steps
        << ")\n"
        << "local SLD response hierarchy test: "
        << (response_hierarchy_ok ? "PASS" : "FAIL")
        << " (depth=" << response_hierarchy.constructed_depth
        << ", projection="
        << static_cast<double>(
               response_hierarchy.final_projection_energy)
        << ", Gram error="
        << static_cast<double>(response_hierarchy.maximum_gram_error)
        << ")\n"
        << "local SLD response family test: "
        << (response_family_ok ? "PASS" : "FAIL")
        << " (rows=" << response_family.rows.size()
        << ", coefficient difference="
        << static_cast<double>(
               response_family.maximum_coefficient_l2_difference)
        << ")\n"
        << "local SLD trajectory evaluator test: "
        << (trajectory_evaluation_ok ? "PASS" : "FAIL")
        << " (maximum="
        << static_cast<double>(
               trajectory_evaluation.maximum.terminal_ratio)
        << ", peak step=" << trajectory_evaluation.maximum.steps
        << ", time refinement="
        << static_cast<double>(
               trajectory_evaluation.time_step_relative_error)
        << ")\n"
        << "local SLD signature block test: "
        << (signature_block_ok ? "PASS" : "FAIL")
        << " (dominant="
        << static_cast<double>(
               signature_block.dominant_closed_sld_ratio)
        << ", cross="
        << static_cast<double>(signature_block.cross_sld_ratio)
        << ", reconstruction error="
        << static_cast<double>(
               signature_block.quotient_reconstruction_error)
        << ")\n"
        << "static adversary test: " << (adversary_ok ? "PASS" : "FAIL")
        << " (Q=" << static_cast<double>(adversary.objective.energy_level_quantity)
        << ")\n"
        << "dynamic adversary ensemble test: "
        << (dynamic_class_ok ? "PASS" : "FAIL")
        << " (restarts=" << dynamic_class_result.restart_objectives.size()
        << ", winner=" << dynamic_class_result.winning_restart << ")\n"
        << "Q directional derivative test: "
        << (q_derivative_ok ? "PASS" : "FAIL")
        << " (refinement error="
        << static_cast<double>(q_derivative.relative_refinement_error) << ")\n"
        << "Galerkin RK4 test: " << (evolution_ok ? "PASS" : "FAIL")
        << " (energy residual="
        << static_cast<double>(evolution.energy_balance_residual) << ")\n";
    return rational_ok && scaling_ok && concentration_ok && strong_l4_ok &&
           shifted_density_ok &&
           dyadic_tail_scaling_ok && dyadic_shell_bounds_ok &&
           periodic_shell_geometry_ok && periodic_tail_bound_ok &&
           far_tail_closure_ok &&
           transition_block_scaling_ok &&
           moving_gap_controller_ok &&
           triad_ok && helical_ok && helical_gap_ok && local_symmetry_ok &&
           orthogonal_geometry_ok && local_signature_geometry_ok &&
           local_signature_objective_ok && pure_helical_ok && fft_ok &&
           helical_sector_objective_ok && helical_adversary_ok &&
           helical_trajectory_adjoint_ok &&
           helical_trajectory_adversary_ok && fft_adjoint_ok && adjoint_ok &&
           q_gradient_ok &&
           trajectory_gradient_ok && q_gain_gradient_ok &&
           q_increase_gradient_ok && q_increase_constraints_ok &&
           critical_integral_gradient_ok &&
           partition_integral_gradients_ok && gradient_search_ok &&
           closure_gradient_search_ok &&
           signed_closure_gradient_search_ok &&
           local_sld_gradient_search_ok &&
           cyclic_ansatz_ok && cyclic_trajectory_ansatz_ok &&
           cyclic_krylov_ansatz_ok && signature_block_ok &&
           response_hierarchy_ok && response_family_ok &&
           trajectory_evaluation_ok &&
           adversary_ok && dynamic_class_ok && q_derivative_ok && evolution_ok;
}

int run_adversary(const AdversaryOptions& options, std::ostream& out) {
    active_galerkin.configure(options.backend, 1);  // restart-level parallelism first
    if (!options.state_directory.empty()) {
        std::filesystem::create_directories(
            std::filesystem::path(options.state_directory) / "static");
        std::filesystem::create_directories(
            std::filesystem::path(options.state_directory) / "dynamic");
    }
    std::vector<AdversaryResult> results;
    std::vector<DynamicAdversaryResult> dynamic_results;
    SpectralState replayed_dynamic_warm_state;
    if (!options.dynamic_warm_state.empty()) {
        replayed_dynamic_warm_state =
            SpectralStateReader::read_tsv(options.dynamic_warm_state);
    }
    const LemmaAdversary adversary(options.threads);
    const DynamicAdversaryEnsemble dynamic_adversary(
        options.backend, adversary.threads());
    DynamicAdversaryOptions dynamic_options;
    dynamic_options.generations = options.dynamic_generations;
    dynamic_options.mutation = static_cast<Real>(options.mutation);
    dynamic_options.viscosity = static_cast<Real>(options.viscosity);
    dynamic_options.final_time = static_cast<Real>(options.evolution_time);
    dynamic_options.time_step = static_cast<Real>(options.time_step);
    dynamic_options.seed = options.seed;
    dynamic_options.objective = options.dynamic_objective;
    dynamic_options.optimizer = options.dynamic_optimizer;
    dynamic_options.gradient_method = options.gradient_method;
    dynamic_options.sobolev_order = options.sobolev_order;
    dynamic_options.sobolev_cap = static_cast<Real>(options.sobolev_cap);
    dynamic_options.critical_density_shift =
        static_cast<Real>(options.critical_density_shift);
    dynamic_options.minimum_dyadic_gap = options.minimum_dyadic_gap;
    for (const int cutoff : options.cutoffs) {
        SpectralState warm_start;
        const SpectralState* warm_start_pointer = nullptr;
        if (!results.empty()) {
            warm_start = results.back().state;
            warm_start_pointer = &warm_start;
        }
        AdversaryResult result = optimize_static_depletion_parallel(
            cutoff, options.restarts, options.generations,
            static_cast<Real>(options.mutation), options.seed, warm_start_pointer,
            adversary);
        if (!options.state_prefix.empty()) {
            write_spectral_state(options.state_prefix + "-K" +
                                     std::to_string(cutoff) + ".tsv",
                                 result);
        }
        if (!options.state_directory.empty()) {
            write_spectral_state(
                (std::filesystem::path(options.state_directory) / "static" /
                 ("K" + std::to_string(cutoff) + ".tsv"))
                    .string(),
                result);
        }
        const SpectralState* dynamic_warm_start = nullptr;
        if (options.dynamic_replay_each_cutoff &&
            !replayed_dynamic_warm_state.waves.empty()) {
            dynamic_warm_start = &replayed_dynamic_warm_state;
        } else if (!dynamic_results.empty()) {
            dynamic_warm_start = &dynamic_results.back().state;
        } else if (!replayed_dynamic_warm_state.waves.empty()) {
            dynamic_warm_start = &replayed_dynamic_warm_state;
        }
        DynamicAdversaryResult dynamic = dynamic_adversary.optimize(
            result.state, dynamic_warm_start, dynamic_options,
            options.dynamic_restarts);
        if (!options.state_prefix.empty()) {
            AdversaryResult dynamic_state;
            dynamic_state.cutoff = cutoff;
            dynamic_state.modes = static_cast<int>(dynamic.state.waves.size());
            dynamic_state.state = dynamic.state;
            dynamic_state.objective = dynamic.initial_objective;
            write_spectral_state(options.state_prefix + "-dynamic-K" +
                                     std::to_string(cutoff) + ".tsv",
                                     dynamic_state);
        }
        if (!options.state_directory.empty()) {
            AdversaryResult dynamic_state;
            dynamic_state.cutoff = cutoff;
            dynamic_state.modes =
                static_cast<int>(dynamic.state.waves.size());
            dynamic_state.state = dynamic.state;
            dynamic_state.objective = dynamic.initial_objective;
            write_spectral_state(
                (std::filesystem::path(options.state_directory) / "dynamic" /
                 ("K" + std::to_string(cutoff) + ".tsv"))
                    .string(),
                dynamic_state);
        }
        dynamic_results.push_back(std::move(dynamic));
        results.push_back(std::move(result));
    }

    Real q_growth_ratio = 1.0L;
    Real q_log_slope = 0.0L;
    if (results.size() >= 2) {
        const auto& low = results.front();
        const auto& high = results.back();
        if (low.objective.energy_level_quantity > 0.0L &&
            high.objective.energy_level_quantity > 0.0L &&
            low.cutoff != high.cutoff) {
            q_growth_ratio = high.objective.energy_level_quantity /
                             low.objective.energy_level_quantity;
            q_log_slope =
                std::log(q_growth_ratio) /
                std::log(static_cast<Real>(high.cutoff) /
                         static_cast<Real>(low.cutoff));
        }
    }
    bool embedding_monotonicity = true;
    for (std::size_t index = 1; index < results.size(); ++index) {
        const Real previous = results[index - 1].objective.energy_level_quantity;
        const Real current = results[index].objective.energy_level_quantity;
        embedding_monotonicity = embedding_monotonicity &&
                                 current + 1e-18L * std::max(1.0L, previous) >= previous;
    }
    AdversaryReport report;
    report.workers = adversary.threads();
    report.backend = options.backend;
    report.dynamic_objective = options.dynamic_objective;
    report.dynamic_optimizer = options.dynamic_optimizer;
    report.gradient_method = options.gradient_method;
    report.minimum_dyadic_gap = options.minimum_dyadic_gap;
    report.sobolev_order = options.sobolev_order;
    report.sobolev_cap = static_cast<Real>(options.sobolev_cap);
    report.critical_density_shift =
        static_cast<Real>(options.critical_density_shift);
    report.restarts = options.restarts;
    report.dynamic_restarts = options.dynamic_restarts;
    report.dynamic_replay_each_cutoff =
        options.dynamic_replay_each_cutoff;
    report.generations = options.generations;
    report.dynamic_generations = options.dynamic_generations;
    report.mutation = static_cast<Real>(options.mutation);
    report.seed = options.seed;
    report.viscosity = static_cast<Real>(options.viscosity);
    report.time = static_cast<Real>(options.evolution_time);
    report.requested_dt = static_cast<Real>(options.time_step);
    report.q_growth_ratio = q_growth_ratio;
    report.q_cutoff_log_slope = q_log_slope;
    report.embedding_monotonicity = embedding_monotonicity;
    report.rows.reserve(results.size());
    for (std::size_t index = 0; index < results.size(); ++index) {
        const auto& result = results[index];
        const auto& dynamic = dynamic_results[index];
        const auto& evolution = dynamic.refined_evolution;
        const Real vortex_partition_denominator =
            evolution.integral_absolute_local_vortex +
            evolution.integral_absolute_nonlocal_vortex;
        const Real nonlocal_vortex_fraction =
            vortex_partition_denominator > 0.0L
                ? evolution.integral_absolute_nonlocal_vortex /
                      vortex_partition_denominator
                : 0.0L;
        const Real strong_l4_envelope =
            evolution.maximum_energy_level_quantity * evolution.initial_energy /
            (2.0L * static_cast<Real>(options.viscosity));
        const Real envelope_utilization = strong_l4_envelope > 0.0L
            ? evolution.integral_critical / strong_l4_envelope
            : 0.0L;
        AdversaryReportRow row;
        row.cutoff = result.cutoff;
        row.modes = result.modes;
        row.evaluations = result.evaluations;
        row.accepted_mutations = result.accepted_mutations;
        row.energy = result.objective.energy;
        row.enstrophy = result.objective.enstrophy;
        row.palinstrophy = result.objective.palinstrophy;
        row.vortex_stretching = result.objective.vortex_stretching;
        row.depletion = result.objective.depletion;
        row.q = result.objective.energy_level_quantity;
        row.critical_integrand = result.objective.critical_integrand;
        row.dynamic_steps = evolution.steps;
        row.dynamic_integral = evolution.integral_critical;
        row.dynamic_coarse_integral = dynamic.evolution.integral_critical;
        row.dynamic_local_integral = evolution.integral_local_critical;
        row.dynamic_nonlocal_integral = evolution.integral_nonlocal_critical;
        row.dynamic_near_nonlocal_integral =
            evolution.integral_near_nonlocal_critical;
        row.dynamic_far_nonlocal_integral =
            evolution.integral_far_nonlocal_critical;
        row.dynamic_selected_gap_tail_integral =
            evolution.integral_selected_gap_tail_critical;
        row.dynamic_dt_relative_error = dynamic.time_step_relative_error;
        row.dynamic_search_initial_objective =
            dynamic.search_initial_objective;
        row.dynamic_search_final_objective =
            dynamic.search_final_objective;
        row.dynamic_initial_local_critical_density =
            evolution.initial_local_critical_integrand;
        row.dynamic_final_local_critical_density =
            evolution.final_local_critical_integrand;
        row.dynamic_initial_enstrophy = evolution.initial_enstrophy;
        row.dynamic_initial_palinstrophy = evolution.initial_palinstrophy;
        row.dynamic_initial_ep_shift =
            evolution.initial_energy * evolution.initial_palinstrophy;
        const Real initial_frequency = std::sqrt(
            evolution.initial_enstrophy /
            std::max(1e-30L, evolution.initial_energy));
        const Real local_growth_normalization = options.evolution_time *
            initial_frequency * evolution.initial_enstrophy;
        if (evolution.initial_local_critical_integrand > 1e-30L &&
            evolution.final_local_critical_integrand > 1e-30L) {
            row.dynamic_local_critical_log_gain = std::log(
                evolution.final_local_critical_integrand /
                evolution.initial_local_critical_integrand);
            if (local_growth_normalization > 1e-30L) {
                row.dynamic_local_log_gain_rate_ratio =
                    row.dynamic_local_critical_log_gain /
                    local_growth_normalization;
            }
        }
        const Real ep_shifted_initial =
            evolution.initial_local_critical_integrand +
            row.dynamic_initial_ep_shift;
        const Real ep_shifted_final =
            evolution.final_local_critical_integrand +
            row.dynamic_initial_ep_shift;
        if (ep_shifted_initial > 1e-30L && ep_shifted_final > 1e-30L) {
            row.dynamic_ep_shifted_local_log_gain =
                std::log1p(
                    (evolution.final_local_critical_integrand -
                     evolution.initial_local_critical_integrand) /
                    ep_shifted_initial);
            if (local_growth_normalization > 1e-30L) {
                row.dynamic_ep_shifted_log_gain_rate_ratio =
                    row.dynamic_ep_shifted_local_log_gain /
                    local_growth_normalization;
            }
        }
        row.dynamic_maximum_q = evolution.maximum_energy_level_quantity;
        row.dynamic_initial_q = evolution.initial_energy_level_quantity;
        row.dynamic_final_q = evolution.final_energy_level_quantity;
        row.dynamic_log_q_gain =
            evolution.initial_energy_level_quantity > 1e-30L &&
                    evolution.final_energy_level_quantity > 1e-30L
                ? std::log(evolution.final_energy_level_quantity /
                           evolution.initial_energy_level_quantity)
                : -std::numeric_limits<Real>::infinity();
        row.dynamic_maximum_local_q =
            evolution.maximum_local_energy_level_quantity;
        row.dynamic_maximum_nonlocal_q =
            evolution.maximum_nonlocal_energy_level_quantity;
        row.dynamic_maximum_near_nonlocal_q =
            evolution.maximum_near_nonlocal_energy_level_quantity;
        row.dynamic_maximum_far_nonlocal_q =
            evolution.maximum_far_nonlocal_energy_level_quantity;
        row.dynamic_maximum_selected_gap_tail_q =
            evolution.maximum_selected_gap_tail_energy_level_quantity;
        row.dynamic_q_log_growth_ratio =
            evolution.maximum_positive_q_log_growth_ratio;
        row.dynamic_q_derivative_error =
            evolution.maximum_q_derivative_refinement_error;
        row.strong_l4_envelope = strong_l4_envelope;
        row.envelope_utilization = envelope_utilization;
        row.dynamic_maximum_enstrophy = evolution.maximum_enstrophy;
        row.dynamic_maximum_vorticity = evolution.maximum_vorticity_linf;
        row.dynamic_maximum_holder_half =
            evolution.maximum_holder_half_coherence;
        row.dynamic_maximum_stretch_alignment =
            evolution.maximum_stretch_alignment;
        row.dynamic_nonlocal_vortex_fraction = nonlocal_vortex_fraction;
        row.dynamic_partition_residual =
            evolution.maximum_vortex_partition_residual;
        row.dynamic_final_energy = evolution.final_energy;
        row.dynamic_energy_balance_residual = evolution.energy_balance_residual;
        row.dynamic_integral_absolute_local_vortex =
            evolution.integral_absolute_local_vortex;
        row.dynamic_integral_absolute_nonlocal_vortex =
            evolution.integral_absolute_nonlocal_vortex;
        row.dynamic_integral_absolute_total_vortex =
            evolution.integral_absolute_total_vortex;
        row.dynamic_geometry_samples = evolution.geometry_samples;
        row.dynamic_evaluations = dynamic.evaluations;
        row.dynamic_winning_restart = dynamic.winning_restart;
        row.dynamic_restart_objectives.assign(
            dynamic.restart_objectives.begin(),
            dynamic.restart_objectives.end());
        row.dynamic_accepted_mutations = dynamic.accepted_mutations;
        row.dynamic_accepted_gradient_steps =
            dynamic.accepted_gradient_steps;
        row.dynamic_sobolev_value = InitialSobolevConstraint(
            options.sobolev_order,
            static_cast<Real>(options.sobolev_cap))
                                         .value(dynamic.state);
        row.dynamic_gradient_trace.reserve(dynamic.gradient_trace.size());
        for (const GradientIterationRecord& point : dynamic.gradient_trace) {
            row.dynamic_gradient_trace.push_back(AdversaryGradientTracePoint{
                point.iteration,
                point.objective_before,
                point.objective_after,
                point.projected_gradient_norm,
                point.accepted_step,
                point.sobolev_value,
                point.line_search_evaluations,
                point.used_steepest_fallback,
                point.accepted});
        }
        report.rows.push_back(row);
    }
    AdversaryReporter::write_console(report, out);
    if (!options.certificate_path.empty()) {
        std::ofstream certificate(options.certificate_path);
        if (!certificate) {
            throw std::runtime_error("cannot open adversary certificate: " +
                                     options.certificate_path);
        }
        AdversaryReporter::write_json(report, certificate);
        out << "Certificate written to " << options.certificate_path << '\n';
    }
    const bool evolutions_valid = std::all_of(
        dynamic_results.begin(), dynamic_results.end(),
        [](const DynamicAdversaryResult& dynamic) {
            const EvolutionResult& evolution = dynamic.refined_evolution;
            return evolution.finite &&
                   evolution.final_energy <= evolution.initial_energy + 1e-12L &&
                   std::abs(evolution.energy_balance_residual) < 1e-6L &&
                   dynamic.time_step_relative_error < 1e-4L;
        });
    return embedding_monotonicity && evolutions_valid &&
           std::all_of(results.begin(), results.end(), [](const AdversaryResult& result) {
        return std::isfinite(result.objective.energy_level_quantity) &&
               std::abs(result.objective.energy - 1.0L) < 1e-12L;
    }) ? 0 : 2;
}

int run_family(const FamilyOptions& options, std::ostream& out) {
    active_galerkin.configure(options.backend, 1);
    struct FamilyRun {
        std::uint64_t seed = 0;
        int cutoff = 0;
        SpectralState initial;
        StaticObjective objective;
        EvolutionResult coarse;
        EvolutionResult refined;
        Real dt_error = 0.0L;
        Real dt_absolute_error = 0.0L;
        Real projection_residual = 0.0L;
        Real q_enstrophy_envelope = 0.0L;
        Real energy_identity_envelope = 0.0L;
        Real envelope_utilization = 0.0L;
        Real factorization_violation = 0.0L;
    };

    std::vector<FamilyRun> runs;
    const std::size_t cutoffs_per_seed = options.cutoffs.size();
    runs.reserve(cutoffs_per_seed * static_cast<std::size_t>(options.seed_count));
    for (int seed_offset = 0; seed_offset < options.seed_count; ++seed_offset) {
        const std::uint64_t seed =
            options.seed + static_cast<std::uint64_t>(seed_offset);
        const std::size_t family_begin = runs.size();
        for (const int cutoff : options.cutoffs) {
            FamilyRun run;
            run.seed = seed;
            run.cutoff = cutoff;
            run.initial = SpectralStateFactory::analytic(
                cutoff, seed, static_cast<Real>(options.spectral_decay));
            runs.push_back(std::move(run));
        }
        const std::size_t family_end = runs.size();
        const Real maximum_cutoff_energy =
            SpectralStateOps::energy(runs[family_end - 1].initial);
        if (!(maximum_cutoff_energy > 0.0L)) {
            throw std::runtime_error("consistent family has zero energy");
        }
        const Real common_factor = 1.0L / std::sqrt(maximum_cutoff_energy);
        for (std::size_t index = family_begin; index < family_end; ++index) {
            SpectralStateOps::scale(runs[index].initial, common_factor);
        }

        // Every lower state must literally be the restriction of the next.
        for (std::size_t index = family_begin + 1; index < family_end; ++index) {
            const SpectralState& lower = runs[index - 1].initial;
            const SpectralState& upper = runs[index].initial;
            for (std::size_t mode = 0; mode < lower.waves.size(); ++mode) {
                const auto upper_mode = upper.index.find(lower.waves[mode]);
                if (upper_mode == upper.index.end()) {
                    throw std::runtime_error("non-nested family cutoff list");
                }
                for (std::size_t direction = 0; direction < 3; ++direction) {
                    runs[index].projection_residual = std::max(
                        runs[index].projection_residual,
                        std::abs(lower.velocity[mode][direction] -
                                 upper.velocity[upper_mode->second][direction]));
                }
            }
        }
    }

    // Build the shared read-only convolution tables before worker threads start.
    for (const auto& run : runs) {
        static_cast<void>(SpectralStateOps::interactions(run.initial));
    }
    const ProjectiveFamily family(options.threads);
    const bool internal_parallelism =
        options.backend == "fft" ||
        (options.backend == "auto" && options.cutoffs.back() >= 5);
    active_galerkin.set_compute_threads(internal_parallelism ? family.threads() : 1);
    auto process_run = [&](std::size_t index) {
        auto& run = runs[index];
        run.objective = active_trajectory_analyzer.evaluate_static(run.initial);
        run.coarse = active_trajectory_analyzer.evolve(
            run.initial, static_cast<Real>(options.viscosity),
            static_cast<Real>(options.evolution_time),
            static_cast<Real>(options.time_step));
        run.refined = active_trajectory_analyzer.evolve(
            run.initial, static_cast<Real>(options.viscosity),
            static_cast<Real>(options.evolution_time),
            0.5L * static_cast<Real>(options.time_step), true);
        run.dt_absolute_error =
            std::abs(run.refined.integral_critical - run.coarse.integral_critical);
        run.dt_error = run.dt_absolute_error /
            std::max(1e-30L, std::abs(run.refined.integral_critical));
        run.q_enstrophy_envelope =
            run.refined.maximum_energy_level_quantity *
            run.refined.integral_enstrophy;
        run.energy_identity_envelope =
            run.refined.maximum_energy_level_quantity *
            run.refined.initial_energy /
            (2.0L * static_cast<Real>(options.viscosity));
        run.envelope_utilization = run.energy_identity_envelope > 0.0L
            ? run.refined.integral_critical / run.energy_identity_envelope
            : 0.0L;
        run.factorization_violation = std::max(
            0.0L, run.refined.integral_critical - run.q_enstrophy_envelope);
    };
    family.run_cutoffs(runs.size(), internal_parallelism, process_run);
    active_galerkin.set_compute_threads(1);

    std::vector<FamilySummaryRow> summaries;
    summaries.reserve(static_cast<std::size_t>(options.seed_count));
    for (int seed_offset = 0; seed_offset < options.seed_count; ++seed_offset) {
        const std::size_t begin =
            static_cast<std::size_t>(seed_offset) * cutoffs_per_seed;
        const std::size_t end = begin + cutoffs_per_seed;
        FamilySummaryRow summary;
        summary.seed = runs[begin].seed;
        if (cutoffs_per_seed >= 2) {
            summary.last_increment = runs[end - 1].refined.integral_critical -
                                     runs[end - 2].refined.integral_critical;
            summary.last_relative_increment =
                std::abs(summary.last_increment) /
                std::max(1e-30L,
                         std::abs(runs[end - 1].refined.integral_critical));
            const Real low_q =
                runs[begin].refined.maximum_energy_level_quantity;
            const Real high_q =
                runs[end - 1].refined.maximum_energy_level_quantity;
            if (low_q > 0.0L && high_q > 0.0L &&
                runs[begin].cutoff != runs[end - 1].cutoff) {
                summary.endpoint_q_growth_ratio = high_q / low_q;
                summary.endpoint_q_log_slope =
                    std::log(summary.endpoint_q_growth_ratio) /
                    std::log(static_cast<Real>(runs[end - 1].cutoff) /
                             static_cast<Real>(runs[begin].cutoff));
            }
            Real previous_running_max = 0.0L;
            for (std::size_t index = begin; index + 1 < end; ++index) {
                previous_running_max = std::max(
                    previous_running_max,
                    runs[index].refined.maximum_energy_level_quantity);
            }
            summary.maximum_q = std::max(previous_running_max, high_q);
            if (previous_running_max > 0.0L && high_q > previous_running_max) {
                summary.tail_record_growth_ratio = high_q / previous_running_max;
                summary.tail_record_log_slope =
                    std::log(summary.tail_record_growth_ratio) /
                    std::log(static_cast<Real>(runs[end - 1].cutoff) /
                             static_cast<Real>(runs[end - 2].cutoff));
            }
        } else {
            summary.maximum_q =
                runs[begin].refined.maximum_energy_level_quantity;
        }
        summaries.push_back(summary);
    }
    const auto worst_q = std::max_element(
        summaries.begin(), summaries.end(),
        [](const FamilySummaryRow& left, const FamilySummaryRow& right) {
            return left.tail_record_log_slope < right.tail_record_log_slope;
        });
    const auto worst_increment = std::max_element(
        summaries.begin(), summaries.end(),
        [](const FamilySummaryRow& left, const FamilySummaryRow& right) {
            return left.last_relative_increment < right.last_relative_increment;
        });

    FamilyReport report;
    report.initial_seed = options.seed;
    report.seed_count = options.seed_count;
    report.spectral_decay = static_cast<Real>(options.spectral_decay);
    report.viscosity = static_cast<Real>(options.viscosity);
    report.time = static_cast<Real>(options.evolution_time);
    report.threads = family.threads();
    report.backend = options.backend;
    report.summaries = summaries;
    report.worst_tail_seed = worst_q->seed;
    report.worst_tail_growth_ratio = worst_q->tail_record_growth_ratio;
    report.worst_tail_log_slope = worst_q->tail_record_log_slope;
    report.worst_last_relative_increment =
        worst_increment->last_relative_increment;
    report.runs.reserve(runs.size());
    for (const auto& run : runs) {
        FamilyReportRow row;
        row.seed = run.seed;
        row.cutoff = run.cutoff;
        row.modes = run.initial.waves.size();
        row.initial_energy = run.objective.energy;
        row.initial_enstrophy = run.objective.enstrophy;
        row.integral_critical = run.refined.integral_critical;
        row.maximum_q = run.refined.maximum_energy_level_quantity;
        row.maximum_local_q = run.refined.maximum_local_energy_level_quantity;
        row.maximum_nonlocal_q = run.refined.maximum_nonlocal_energy_level_quantity;
        row.maximum_positive_q_log_growth_ratio =
            run.refined.maximum_positive_q_log_growth_ratio;
        row.q_derivative_refinement_error =
            run.refined.maximum_q_derivative_refinement_error;
        row.q_derivative_samples = run.refined.q_derivative_samples;
        row.q_enstrophy_envelope = run.q_enstrophy_envelope;
        row.energy_identity_envelope = run.energy_identity_envelope;
        row.envelope_utilization = run.envelope_utilization;
        row.factorization_violation = run.factorization_violation;
        row.dt_absolute_error = run.dt_absolute_error;
        row.dt_relative_error = run.dt_error;
        row.maximum_enstrophy = run.refined.maximum_enstrophy;
        row.maximum_vorticity = run.refined.maximum_vorticity_linf;
        row.maximum_holder_half = run.refined.maximum_holder_half_coherence;
        row.local_integral = run.refined.integral_local_critical;
        row.nonlocal_integral = run.refined.integral_nonlocal_critical;
        row.projection_residual = run.projection_residual;
        report.runs.push_back(row);
    }
    FamilyReporter::write_console(report, out);
    if (!options.certificate_path.empty()) {
        std::ofstream certificate(options.certificate_path);
        if (!certificate) {
            throw std::runtime_error("cannot open family certificate: " +
                                     options.certificate_path);
        }
        FamilyReporter::write_json(report, certificate);
        out << "Certificate written to " << options.certificate_path << '\n';
    }

    const bool valid = std::all_of(runs.begin(), runs.end(), [](const FamilyRun& run) {
        return run.refined.finite && run.projection_residual < 1e-18L &&
               (run.dt_error < 1e-4L || run.dt_absolute_error < 1e-12L) &&
               run.factorization_violation < 1e-15L &&
               std::abs(run.refined.energy_balance_residual) < 1e-6L;
    });
    return valid ? 0 : 2;
}

}  // namespace lemma
