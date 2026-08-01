#include "lemma_engine.hpp"
#include "adversary_reporter.hpp"
#include "family_reporter.hpp"
#include "gradient_adversary.hpp"
#include "proof_scaling.hpp"
#include "parallel_executor.hpp"
#include "lemma_adversary.hpp"
#include "lemma_reporter.hpp"
#include "projective_family.hpp"
#include "spectral_adjoint.hpp"
#include "spectral_dynamics.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_objective.hpp"
#include "spectral_state.hpp"
#include "state_analysis.hpp"
#include "trajectory_analyzer.hpp"
#include "triad_verifier.hpp"

#include <algorithm>
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

struct DynamicAdversaryResult {
    SpectralState state;
    StaticObjective initial_objective;
    EvolutionResult evolution;
    EvolutionResult refined_evolution;
    Real time_step_relative_error = 0.0L;
    Real search_initial_objective = 0.0L;
    Real search_final_objective = 0.0L;
    std::vector<GradientIterationRecord> gradient_trace;
    int accepted_mutations = 0;
    int accepted_gradient_steps = 0;
    int evaluations = 0;
};

Real dynamic_objective_value(const EvolutionResult& evolution,
                             const std::string& objective) {
    if (objective == "critical-integral") {
        return evolution.integral_critical;
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
            return -std::numeric_limits<Real>::infinity();
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

DynamicAdversaryResult optimize_dynamic(
    const SpectralState& primary_start, const SpectralState* secondary_start,
    int generations, Real mutation, Real viscosity, Real final_time, Real dt,
    std::uint64_t seed, const std::string& objective,
    const std::string& optimizer, const std::string& gradient_method,
    int sobolev_order, Real sobolev_cap) {
    if (generations < 0) {
        throw std::invalid_argument("--dynamic-generations cannot be negative");
    }
    std::mt19937_64 generator(seed ^ 0xd1b54a32d192ed03ULL ^
                              static_cast<std::uint64_t>(SpectralStateOps::cutoff(primary_start)));
    DynamicAdversaryResult result;
    const InitialSobolevConstraint sobolev(sobolev_order, sobolev_cap);
    result.state = primary_start;
    result.initial_objective = active_trajectory_analyzer.evaluate_static(result.state);
    result.evolution = active_trajectory_analyzer.evolve(result.state, viscosity, final_time, dt);
    ++result.evaluations;
    bool result_admissible = sobolev.admissible(result.state);

    if (secondary_start != nullptr) {
        const int target_cutoff = SpectralStateOps::cutoff(primary_start);
        const int source_cutoff = SpectralStateOps::cutoff(*secondary_start);
        SpectralState secondary = source_cutoff <= target_cutoff
            ? SpectralStateFactory::lift(
                  *secondary_start, target_cutoff, generator)
            : SpectralStateFactory::project(
                  *secondary_start, target_cutoff);
        const EvolutionResult secondary_evolution =
            active_trajectory_analyzer.evolve(secondary, viscosity, final_time, dt);
        ++result.evaluations;
        const bool secondary_admissible = sobolev.admissible(secondary);
        if (secondary_admissible &&
            (!result_admissible ||
             dynamic_objective_value(secondary_evolution, objective) >
                 dynamic_objective_value(result.evolution, objective))) {
            result.state = std::move(secondary);
            result.initial_objective = active_trajectory_analyzer.evaluate_static(result.state);
            result.evolution = secondary_evolution;
            result_admissible = true;
        }
    }
    if (!result_admissible) {
        throw std::invalid_argument(
            "no dynamic start satisfies the configured Sobolev cap");
    }
    result.search_initial_objective =
        dynamic_objective_value(result.evolution, objective);

    const bool use_mutations = optimizer == "mutate" || optimizer == "hybrid";
    const bool use_gradient = optimizer == "gradient" || optimizer == "hybrid";
    for (int generation = 0; use_mutations && generation < generations;
         ++generation) {
        const Real progress = generations > 1
                                  ? static_cast<Real>(generation) /
                                        static_cast<Real>(generations - 1)
                                  : 0.0L;
        const Real radius = mutation * (0.5L - 0.4L * progress);
        SpectralState candidate =
            SpectralStateFactory::mutate(
                result.state, radius, generator, generation % 4 != 0);
        if (!sobolev.admissible(candidate)) {
            continue;
        }
        const EvolutionResult evolution =
            active_trajectory_analyzer.evolve(candidate, viscosity, final_time, dt);
        ++result.evaluations;
        if (evolution.finite &&
            dynamic_objective_value(evolution, objective) >
                dynamic_objective_value(result.evolution, objective)) {
            result.state = std::move(candidate);
            result.initial_objective = active_trajectory_analyzer.evaluate_static(result.state);
            result.evolution = evolution;
            ++result.accepted_mutations;
        }
    }
    if (use_gradient && generations > 0) {
        const int trajectory_steps = std::max(
            1, static_cast<int>(std::ceil(final_time / dt)));
        GradientSearchOptions gradient_options;
        gradient_options.iterations = generations;
        gradient_options.line_search_steps = 16;
        gradient_options.trajectory_steps = trajectory_steps;
        gradient_options.viscosity = viscosity;
        gradient_options.time_step =
            final_time / static_cast<Real>(trajectory_steps);
        gradient_options.initial_step = mutation;
        gradient_options.objective = objective;
        gradient_options.method = gradient_method;
        gradient_options.sobolev_order = sobolev_order;
        gradient_options.sobolev_cap = sobolev_cap;
        const GradientSearchResult gradient =
            active_gradient_adversary.maximize_q(
                result.state, gradient_options);
        result.state = gradient.state;
        result.initial_objective = active_trajectory_analyzer.evaluate_static(result.state);
        result.evolution =
            active_trajectory_analyzer.evolve(result.state, viscosity, final_time, dt);
        result.evaluations += gradient.trajectory_evaluations + 1;
        result.accepted_gradient_steps = gradient.accepted_steps;
        result.gradient_trace = gradient.trace;
    }
    result.search_final_objective =
        dynamic_objective_value(result.evolution, objective);
    result.refined_evolution =
        active_trajectory_analyzer.evolve(result.state, viscosity, final_time, 0.5L * dt, true);
    result.time_step_relative_error =
        std::abs(result.refined_evolution.integral_critical -
                 result.evolution.integral_critical) /
        std::max(1e-30L, std::abs(result.refined_evolution.integral_critical));
    return result;
}

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
    std::ofstream state_file(path);
    if (!state_file) {
        throw std::runtime_error("cannot open adversarial state file: " + path);
    }
    state_file << "# cutoff=" << result.cutoff << " energy=" << std::setprecision(20)
               << static_cast<double>(result.objective.energy)
               << " Q=D^4*Z=" << static_cast<double>(result.objective.energy_level_quantity)
               << '\n'
               << "kx\tky\tkz\tux_re\tux_im\tuy_re\tuy_im\tuz_re\tuz_im\n";
    for (std::size_t index = 0; index < result.state.waves.size(); ++index) {
        const WaveVector wave = result.state.waves[index];
        state_file << wave.x << '\t' << wave.y << '\t' << wave.z;
        for (const Complex component : result.state.velocity[index]) {
            state_file << '\t' << static_cast<double>(component.real()) << '\t'
                       << static_cast<double>(component.imag());
        }
        state_file << '\n';
    }
}

}  // namespace

int run(const Options& options, std::ostream& out) {
    const ScalingCertificate scaling =
        ScalingAnalyzer::analyze_monomials(options.exponent_denominator);
    const ConcentrationScaling concentration =
        ScalingAnalyzer::analyze_concentration();
    const StrongL4Reduction strong_l4 =
        ScalingAnalyzer::analyze_strong_l4_reduction();
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
    const TriadCertificate triads = TriadVerifier::analyze(2, 2, 7);
    const bool triad_ok = triads.maximum_normalized_energy_residual < 1e-15L &&
                          triads.maximum_divergence_residual < 1e-15L &&
                          triads.maximum_detailed_triad_residual < 1e-15L &&
                          triads.maximum_flux_partition_residual < 1e-15L &&
                          triads.nonzero_vortex_stretching_seen;
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
    const AdversaryResult adversary =
        optimize_static_depletion(1, 1, 2, 0.1L, 11);
    const bool adversary_ok = adversary.modes == 26 &&
                              std::abs(adversary.objective.energy - 1.0L) < 1e-15L &&
                              std::isfinite(adversary.objective.energy_level_quantity) &&
                              adversary.objective.energy_level_quantity >= 0.0L;
    const EvolutionResult evolution =
        active_trajectory_analyzer.evolve(adversary.state, 0.1L, 0.002L, 0.001L);
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
        << "spectral triad test: " << (triad_ok ? "PASS" : "FAIL")
        << " (energy residual="
        << static_cast<double>(triads.maximum_normalized_energy_residual) << ")\n"
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
        << "projected gradient adversary test: "
        << (gradient_search_ok ? "PASS" : "FAIL")
        << " (Q " << static_cast<double>(gradient_search.initial_objective)
        << " -> " << static_cast<double>(gradient_search.objective)
        << ", accepted=" << gradient_search.accepted_steps
        << ", energy error=" << static_cast<double>(gradient_energy_error)
        << ", constraint error="
        << static_cast<double>(gradient_constraint_error)
        << ")\n"
        << "static adversary test: " << (adversary_ok ? "PASS" : "FAIL")
        << " (Q=" << static_cast<double>(adversary.objective.energy_level_quantity)
        << ")\n"
        << "Q directional derivative test: "
        << (q_derivative_ok ? "PASS" : "FAIL")
        << " (refinement error="
        << static_cast<double>(q_derivative.relative_refinement_error) << ")\n"
        << "Galerkin RK4 test: " << (evolution_ok ? "PASS" : "FAIL")
        << " (energy residual="
        << static_cast<double>(evolution.energy_balance_residual) << ")\n";
    return rational_ok && scaling_ok && concentration_ok && strong_l4_ok &&
           triad_ok && fft_ok && fft_adjoint_ok && adjoint_ok && q_gradient_ok &&
           trajectory_gradient_ok && q_gain_gradient_ok &&
           q_increase_gradient_ok && q_increase_constraints_ok &&
           critical_integral_gradient_ok && gradient_search_ok &&
           adversary_ok && q_derivative_ok && evolution_ok;
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
        if (!dynamic_results.empty()) {
            dynamic_warm_start = &dynamic_results.back().state;
        } else if (!replayed_dynamic_warm_state.waves.empty()) {
            dynamic_warm_start = &replayed_dynamic_warm_state;
        }
        active_galerkin.set_compute_threads(adversary.threads());
        DynamicAdversaryResult dynamic = optimize_dynamic(
            result.state, dynamic_warm_start, options.dynamic_generations,
            static_cast<Real>(options.mutation),
            static_cast<Real>(options.viscosity),
            static_cast<Real>(options.evolution_time),
            static_cast<Real>(options.time_step), options.seed,
            options.dynamic_objective, options.dynamic_optimizer,
            options.gradient_method,
            options.sobolev_order,
            static_cast<Real>(options.sobolev_cap));
        active_galerkin.set_compute_threads(1);
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
    report.sobolev_order = options.sobolev_order;
    report.sobolev_cap = static_cast<Real>(options.sobolev_cap);
    report.restarts = options.restarts;
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
        row.dynamic_dt_relative_error = dynamic.time_step_relative_error;
        row.dynamic_search_initial_objective =
            dynamic.search_initial_objective;
        row.dynamic_search_final_objective =
            dynamic.search_final_objective;
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
