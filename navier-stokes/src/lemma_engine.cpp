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

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace lemma {
namespace {

using Integer = SpectralInteger;
using Real = SpectralReal;
using Complex = SpectralComplex;
SpectralGalerkin active_galerkin;
SpectralDynamics active_dynamics(active_galerkin);
SpectralObjective active_objective(active_dynamics);
SpectralAdjoint active_adjoint(active_dynamics, active_objective);
GradientAdversary active_gradient_adversary(
    active_dynamics, active_objective, active_adjoint);

struct VortexPartition {
    Real local = 0.0L;
    Real nonlocal = 0.0L;
    Real absolute_local_pairs = 0.0L;
    Real absolute_nonlocal_pairs = 0.0L;
};

Real critical_integrand_from_stretching(Real stretching, Real enstrophy,
                                        Real palinstrophy) {
    const Real denominator = enstrophy * std::pow(palinstrophy, 3.0L);
    if (!(denominator > 0.0L)) {
        return 0.0L;
    }
    return std::pow(std::abs(stretching), 4.0L) / denominator;
}

Real energy_level_quantity_from_stretching(Real stretching, Real enstrophy,
                                            Real palinstrophy) {
    const Real denominator = std::pow(enstrophy, 2.0L) *
                             std::pow(palinstrophy, 3.0L);
    if (!(denominator > 0.0L)) {
        return 0.0L;
    }
    return std::pow(std::abs(stretching), 4.0L) / denominator;
}

VortexPartition evaluate_vortex_partition(const SpectralState& state) {
    VortexPartition result;
    const Complex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, target_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[target_index];
        const Complex coefficient =
            imaginary_unit * wave_dot(q, state.velocity[p_index]);
        ComplexVector pair{};
        for (std::size_t direction = 0; direction < 3; ++direction) {
            pair[direction] = coefficient * state.velocity[q_index][direction];
        }
        const Real enstrophy_transfer =
            static_cast<Real>(norm_squared(k)) *
            std::real(dot_hermitian(state.velocity[target_index], pair));
        const Integer smallest = std::min(
            {norm_squared(k), norm_squared(p), norm_squared(q)});
        const Integer largest = std::max(
            {norm_squared(k), norm_squared(p), norm_squared(q)});
        if (largest <= 4 * smallest) {
            result.local += enstrophy_transfer;
            result.absolute_local_pairs += std::abs(enstrophy_transfer);
        } else {
            result.nonlocal += enstrophy_transfer;
            result.absolute_nonlocal_pairs += std::abs(enstrophy_transfer);
        }
    }
    return result;
}

using RealVector = std::array<Real, 3>;

struct GeometryDiagnostic {
    Real maximum_vorticity = 0.0L;
    Real rms_vorticity = 0.0L;
    Real maximum_holder_half_coherence = 0.0L;
    Real mean_holder_half_coherence = 0.0L;
    Real positive_stretching_fraction = 0.0L;
    Real maximum_stretch_alignment = 0.0L;
    int high_vorticity_points = 0;
    int coherence_pairs = 0;
};

Real real_vector_norm(const RealVector& vector) {
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] +
                     vector[2] * vector[2]);
}

GeometryDiagnostic evaluate_geometry(const SpectralState& state,
                                     Real high_vorticity_fraction = 0.5L) {
    const int cutoff = SpectralStateOps::cutoff(state);
    const int grid_n = 2 * cutoff + 3;
    const Real spacing = 2.0L * std::acos(-1.0L) / static_cast<Real>(grid_n);
    const std::size_t cells = static_cast<std::size_t>(grid_n) *
                              static_cast<std::size_t>(grid_n) *
                              static_cast<std::size_t>(grid_n);
    std::vector<RealVector> vorticity(cells);
    std::vector<Real> magnitude(cells, 0.0L);
    GeometryDiagnostic result;

    auto cell_index = [grid_n](int x, int y, int z) {
        const auto wrap = [grid_n](int value) {
            value %= grid_n;
            return value < 0 ? value + grid_n : value;
        };
        const auto xx = static_cast<std::size_t>(wrap(x));
        const auto yy = static_cast<std::size_t>(wrap(y));
        const auto zz = static_cast<std::size_t>(wrap(z));
        const auto nn = static_cast<std::size_t>(grid_n);
        return (zz * nn + yy) * nn + xx;
    };

    Real vorticity2_sum = 0.0L;
    Real positive_stretching = 0.0L;
    Real absolute_stretching = 0.0L;
    Real maximum_vorticity = 0.0L;
    Real maximum_stretch_alignment = 0.0L;
    const Complex imaginary_unit{0.0L, 1.0L};
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for collapse(3) num_threads(active_galerkin.compute_threads()) schedule(static) \
    reduction(+ : vorticity2_sum, positive_stretching, absolute_stretching) \
    reduction(max : maximum_vorticity, maximum_stretch_alignment)
#endif
    for (int z_index = 0; z_index < grid_n; ++z_index) {
        for (int y_index = 0; y_index < grid_n; ++y_index) {
            for (int x_index = 0; x_index < grid_n; ++x_index) {
                std::array<std::array<Real, 3>, 3> gradient{};
                const Real x = spacing * static_cast<Real>(x_index);
                const Real y = spacing * static_cast<Real>(y_index);
                const Real z = spacing * static_cast<Real>(z_index);
                for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
                    const WaveVector wave = state.waves[mode];
                    const Real angle = static_cast<Real>(wave.x) * x +
                                       static_cast<Real>(wave.y) * y +
                                       static_cast<Real>(wave.z) * z;
                    const Complex phase{std::cos(angle), std::sin(angle)};
                    const std::array<int, 3> wave_component{wave.x, wave.y, wave.z};
                    for (std::size_t velocity_component = 0;
                         velocity_component < 3; ++velocity_component) {
                        for (std::size_t derivative = 0; derivative < 3; ++derivative) {
                            gradient[velocity_component][derivative] += std::real(
                                imaginary_unit *
                                static_cast<Real>(wave_component[derivative]) *
                                state.velocity[mode][velocity_component] * phase);
                        }
                    }
                }
                const RealVector omega{
                    gradient[2][1] - gradient[1][2],
                    gradient[0][2] - gradient[2][0],
                    gradient[1][0] - gradient[0][1]};
                std::array<std::array<Real, 3>, 3> strain{};
                Real strain2 = 0.0L;
                for (std::size_t row = 0; row < 3; ++row) {
                    for (std::size_t column = 0; column < 3; ++column) {
                        strain[row][column] =
                            0.5L * (gradient[row][column] + gradient[column][row]);
                        strain2 += strain[row][column] * strain[row][column];
                    }
                }
                Real stretch = 0.0L;
                for (std::size_t row = 0; row < 3; ++row) {
                    for (std::size_t column = 0; column < 3; ++column) {
                        stretch += omega[row] * strain[row][column] * omega[column];
                    }
                }
                const auto cell = cell_index(x_index, y_index, z_index);
                vorticity[cell] = omega;
                magnitude[cell] = real_vector_norm(omega);
                maximum_vorticity = std::max(maximum_vorticity, magnitude[cell]);
                vorticity2_sum += magnitude[cell] * magnitude[cell];
                positive_stretching += std::max(0.0L, stretch);
                absolute_stretching += std::abs(stretch);
                if (magnitude[cell] > 0.0L && strain2 > 0.0L) {
                    maximum_stretch_alignment = std::max(
                        maximum_stretch_alignment,
                        std::abs(stretch) /
                            (magnitude[cell] * magnitude[cell] * std::sqrt(strain2)));
                }
            }
        }
    }
    result.maximum_vorticity = maximum_vorticity;
    result.maximum_stretch_alignment = maximum_stretch_alignment;
    result.rms_vorticity = std::sqrt(vorticity2_sum / static_cast<Real>(cells));
    if (absolute_stretching > 0.0L) {
        result.positive_stretching_fraction =
            positive_stretching / absolute_stretching;
    }

    const Real threshold = high_vorticity_fraction * result.maximum_vorticity;
    Real coherence_sum = 0.0L;
    Real maximum_holder_half = 0.0L;
    int high_vorticity_points = 0;
    int coherence_pairs = 0;
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for collapse(3) num_threads(active_galerkin.compute_threads()) schedule(static) \
    reduction(+ : coherence_sum, high_vorticity_points, coherence_pairs) \
    reduction(max : maximum_holder_half)
#endif
    for (int z = 0; z < grid_n; ++z) {
        for (int y = 0; y < grid_n; ++y) {
            for (int x = 0; x < grid_n; ++x) {
                const auto cell = cell_index(x, y, z);
                if (magnitude[cell] < threshold || magnitude[cell] <= 0.0L) {
                    continue;
                }
                ++high_vorticity_points;
                const std::array<std::array<int, 3>, 3> neighbors{{
                    {x + 1, y, z}, {x, y + 1, z}, {x, y, z + 1}}};
                for (const auto neighbor_coordinates : neighbors) {
                    const auto neighbor = cell_index(neighbor_coordinates[0],
                                                     neighbor_coordinates[1],
                                                     neighbor_coordinates[2]);
                    if (magnitude[neighbor] < threshold || magnitude[neighbor] <= 0.0L) {
                        continue;
                    }
                    const RealVector& a = vorticity[cell];
                    const RealVector& b = vorticity[neighbor];
                    const RealVector cross{
                        a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0]};
                    const Real sine = std::min(
                        1.0L, real_vector_norm(cross) /
                                  (magnitude[cell] * magnitude[neighbor]));
                    const Real holder = sine / std::sqrt(spacing);
                    maximum_holder_half = std::max(maximum_holder_half, holder);
                    coherence_sum += holder;
                    ++coherence_pairs;
                }
            }
        }
    }
    result.maximum_holder_half_coherence = maximum_holder_half;
    result.high_vorticity_points = high_vorticity_points;
    result.coherence_pairs = coherence_pairs;
    if (result.coherence_pairs > 0) {
        result.mean_holder_half_coherence =
            coherence_sum / static_cast<Real>(result.coherence_pairs);
    }
    return result;
}

StaticObjective evaluate_static_objective(const SpectralState& state) {
    return active_objective.evaluate(state);
}

struct QDerivativeDiagnostic {
    Real derivative = 0.0L;
    Real log_derivative = 0.0L;
    Real relative_refinement_error = 0.0L;
    bool valid = false;
};

QDerivativeDiagnostic evaluate_q_derivative(const SpectralState& state,
                                            Real viscosity) {
    const StaticObjective center = evaluate_static_objective(state);
    const SpectralIncrement rhs = active_dynamics.rhs(state, viscosity);
    Real rhs_norm2 = 0.0L;
    for (const ComplexVector& value : rhs) {
        rhs_norm2 += std::real(dot_hermitian(value, value));
    }
    QDerivativeDiagnostic result;
    if (!(center.energy > 0.0L) || !(rhs_norm2 > 0.0L) ||
        !(center.energy_level_quantity > 1e-30L)) {
        return result;
    }
    const Real base_step =
        1e-4L * std::sqrt(center.energy / rhs_norm2);
    auto q_at_offset = [&](Real offset) {
        SpectralState shifted = state;
        for (std::size_t mode = 0; mode < shifted.velocity.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                shifted.velocity[mode][component] += offset * rhs[mode][component];
            }
        }
        return evaluate_static_objective(shifted).energy_level_quantity;
    };
    auto central_derivative = [&](Real step) {
        return (q_at_offset(step) - q_at_offset(-step)) / (2.0L * step);
    };
    const Real coarse = central_derivative(base_step);
    const Real refined = central_derivative(0.5L * base_step);
    result.derivative = refined;
    result.log_derivative = refined / center.energy_level_quantity;
    result.relative_refinement_error =
        std::abs(refined - coarse) /
        std::max(1e-30L, std::max(std::abs(refined),
                                 center.energy_level_quantity));
    result.valid = std::isfinite(result.derivative) &&
                   std::isfinite(result.log_derivative) &&
                   std::isfinite(result.relative_refinement_error);
    return result;
}

struct EvolutionResult {
    int steps = 0;
    Real time = 0.0L;
    Real initial_energy = 0.0L;
    Real final_energy = 0.0L;
    Real initial_enstrophy = 0.0L;
    Real final_enstrophy = 0.0L;
    Real initial_energy_level_quantity = 0.0L;
    Real final_energy_level_quantity = 0.0L;
    Real integral_critical = 0.0L;
    Real integral_enstrophy = 0.0L;
    Real maximum_energy_level_quantity = 0.0L;
    Real maximum_critical_integrand = 0.0L;
    Real maximum_enstrophy = 0.0L;
    Real energy_balance_residual = 0.0L;
    Real integral_absolute_local_vortex = 0.0L;
    Real integral_absolute_nonlocal_vortex = 0.0L;
    Real integral_absolute_total_vortex = 0.0L;
    Real integral_local_critical = 0.0L;
    Real integral_nonlocal_critical = 0.0L;
    Real maximum_local_energy_level_quantity = 0.0L;
    Real maximum_nonlocal_energy_level_quantity = 0.0L;
    Real maximum_vortex_partition_residual = 0.0L;
    Real maximum_vorticity_linf = 0.0L;
    Real maximum_holder_half_coherence = 0.0L;
    Real maximum_stretch_alignment = 0.0L;
    Real maximum_positive_q_log_growth_ratio = 0.0L;
    Real maximum_q_derivative_refinement_error = 0.0L;
    int geometry_samples = 0;
    int q_derivative_samples = 0;
    bool finite = true;
};

EvolutionResult evolve_galerkin(SpectralState state, Real viscosity,
                                Real final_time, Real requested_dt,
                                bool collect_vortex_partition = false) {
    if (!(viscosity > 0.0L) || !(final_time > 0.0L) || !(requested_dt > 0.0L)) {
        throw std::invalid_argument("evolution viscosity, time, and dt must be positive");
    }
    EvolutionResult result;
    StaticObjective before = evaluate_static_objective(state);
    result.initial_energy = before.energy;
    result.initial_enstrophy = before.enstrophy;
    result.initial_energy_level_quantity = before.energy_level_quantity;
    result.maximum_energy_level_quantity = before.energy_level_quantity;
    result.maximum_critical_integrand = before.critical_integrand;
    result.maximum_enstrophy = before.enstrophy;
    const Real initial_frequency =
        std::sqrt(before.enstrophy / std::max(1e-30L, before.energy));
    auto sample_q_derivative = [&](const SpectralState& sample_state,
                                   const StaticObjective& sample_objective) {
        const QDerivativeDiagnostic derivative =
            evaluate_q_derivative(sample_state, viscosity);
        if (!derivative.valid) {
            return;
        }
        const Real denominator =
            initial_frequency * sample_objective.enstrophy;
        if (denominator > 0.0L) {
            result.maximum_positive_q_log_growth_ratio = std::max(
                result.maximum_positive_q_log_growth_ratio,
                std::max(0.0L, derivative.log_derivative) / denominator);
        }
        result.maximum_q_derivative_refinement_error = std::max(
            result.maximum_q_derivative_refinement_error,
            derivative.relative_refinement_error);
        ++result.q_derivative_samples;
    };
    VortexPartition partition_before;
    if (collect_vortex_partition) {
        partition_before = evaluate_vortex_partition(state);
        result.maximum_local_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.local, before.enstrophy, before.palinstrophy);
        result.maximum_nonlocal_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.nonlocal, before.enstrophy, before.palinstrophy);
        result.maximum_vortex_partition_residual =
            std::abs(before.vortex_stretching -
                     std::abs(partition_before.local + partition_before.nonlocal));
        const GeometryDiagnostic geometry = evaluate_geometry(state);
        result.maximum_vorticity_linf = geometry.maximum_vorticity;
        result.maximum_holder_half_coherence =
            geometry.maximum_holder_half_coherence;
        result.maximum_stretch_alignment = geometry.maximum_stretch_alignment;
        ++result.geometry_samples;
        sample_q_derivative(state, before);
    }

    const int cutoff = SpectralStateOps::cutoff(state);
    const Real maximum_wave2 = 3.0L * static_cast<Real>(cutoff * cutoff);
    const Real diffusion_dt = 0.5L / (viscosity * maximum_wave2);
    Real time = 0.0L;
    while (time < final_time) {
        const Real dt = std::min({requested_dt, diffusion_dt, final_time - time});
        active_dynamics.rk4_step(state, viscosity, dt);
        const StaticObjective after = evaluate_static_objective(state);
        if (collect_vortex_partition) {
            const VortexPartition partition_after = evaluate_vortex_partition(state);
            const Real local_critical_before = critical_integrand_from_stretching(
                partition_before.local, before.enstrophy, before.palinstrophy);
            const Real local_critical_after = critical_integrand_from_stretching(
                partition_after.local, after.enstrophy, after.palinstrophy);
            const Real nonlocal_critical_before = critical_integrand_from_stretching(
                partition_before.nonlocal, before.enstrophy, before.palinstrophy);
            const Real nonlocal_critical_after = critical_integrand_from_stretching(
                partition_after.nonlocal, after.enstrophy, after.palinstrophy);
            const Real local_q_after = energy_level_quantity_from_stretching(
                partition_after.local, after.enstrophy, after.palinstrophy);
            const Real nonlocal_q_after = energy_level_quantity_from_stretching(
                partition_after.nonlocal, after.enstrophy, after.palinstrophy);
            result.integral_local_critical +=
                0.5L * dt * (local_critical_before + local_critical_after);
            result.integral_nonlocal_critical +=
                0.5L * dt * (nonlocal_critical_before + nonlocal_critical_after);
            result.maximum_local_energy_level_quantity =
                std::max(result.maximum_local_energy_level_quantity, local_q_after);
            result.maximum_nonlocal_energy_level_quantity =
                std::max(result.maximum_nonlocal_energy_level_quantity,
                         nonlocal_q_after);
            result.integral_absolute_local_vortex +=
                0.5L * dt * (std::abs(partition_before.local) +
                             std::abs(partition_after.local));
            result.integral_absolute_nonlocal_vortex +=
                0.5L * dt * (std::abs(partition_before.nonlocal) +
                             std::abs(partition_after.nonlocal));
            result.integral_absolute_total_vortex +=
                0.5L * dt *
                (std::abs(partition_before.local + partition_before.nonlocal) +
                 std::abs(partition_after.local + partition_after.nonlocal));
            result.maximum_vortex_partition_residual =
                std::max(result.maximum_vortex_partition_residual,
                         std::abs(after.vortex_stretching -
                                  std::abs(partition_after.local +
                                           partition_after.nonlocal)));
            partition_before = partition_after;
            if (result.steps % 10 == 0 || time + dt >= final_time) {
                const GeometryDiagnostic geometry = evaluate_geometry(state);
                result.maximum_vorticity_linf =
                    std::max(result.maximum_vorticity_linf,
                             geometry.maximum_vorticity);
                result.maximum_holder_half_coherence =
                    std::max(result.maximum_holder_half_coherence,
                             geometry.maximum_holder_half_coherence);
                result.maximum_stretch_alignment =
                    std::max(result.maximum_stretch_alignment,
                             geometry.maximum_stretch_alignment);
                ++result.geometry_samples;
                sample_q_derivative(state, after);
            }
        }
        result.integral_critical +=
            0.5L * dt * (before.critical_integrand + after.critical_integrand);
        result.integral_enstrophy +=
            0.5L * dt * (before.enstrophy + after.enstrophy);
        result.maximum_energy_level_quantity =
            std::max(result.maximum_energy_level_quantity,
                     after.energy_level_quantity);
        result.maximum_critical_integrand =
            std::max(result.maximum_critical_integrand,
                     after.critical_integrand);
        result.maximum_enstrophy =
            std::max(result.maximum_enstrophy, after.enstrophy);
        result.finite = result.finite && std::isfinite(after.energy) &&
                        std::isfinite(after.enstrophy) &&
                        std::isfinite(after.critical_integrand);
        ++result.steps;
        time += dt;
        before = after;
        if (!result.finite) {
            break;
        }
    }
    result.time = time;
    result.final_energy = before.energy;
    result.final_enstrophy = before.enstrophy;
    result.final_energy_level_quantity = before.energy_level_quantity;
    result.energy_balance_residual = result.final_energy - result.initial_energy +
                                     2.0L * viscosity * result.integral_enstrophy;
    return result;
}

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
    const std::string& optimizer, int sobolev_order, Real sobolev_cap) {
    if (generations < 0) {
        throw std::invalid_argument("--dynamic-generations cannot be negative");
    }
    std::mt19937_64 generator(seed ^ 0xd1b54a32d192ed03ULL ^
                              static_cast<std::uint64_t>(SpectralStateOps::cutoff(primary_start)));
    DynamicAdversaryResult result;
    const InitialSobolevConstraint sobolev(sobolev_order, sobolev_cap);
    result.state = primary_start;
    result.initial_objective = evaluate_static_objective(result.state);
    result.evolution = evolve_galerkin(result.state, viscosity, final_time, dt);
    ++result.evaluations;
    bool result_admissible = sobolev.admissible(result.state);

    if (secondary_start != nullptr) {
        SpectralState secondary =
            SpectralStateFactory::lift(
                *secondary_start, SpectralStateOps::cutoff(primary_start), generator);
        const EvolutionResult secondary_evolution =
            evolve_galerkin(secondary, viscosity, final_time, dt);
        ++result.evaluations;
        const bool secondary_admissible = sobolev.admissible(secondary);
        if (secondary_admissible &&
            (!result_admissible ||
             dynamic_objective_value(secondary_evolution, objective) >
                 dynamic_objective_value(result.evolution, objective))) {
            result.state = std::move(secondary);
            result.initial_objective = evaluate_static_objective(result.state);
            result.evolution = secondary_evolution;
            result_admissible = true;
        }
    }
    if (!result_admissible) {
        throw std::invalid_argument(
            "no dynamic start satisfies the configured Sobolev cap");
    }

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
            evolve_galerkin(candidate, viscosity, final_time, dt);
        ++result.evaluations;
        if (evolution.finite &&
            dynamic_objective_value(evolution, objective) >
                dynamic_objective_value(result.evolution, objective)) {
            result.state = std::move(candidate);
            result.initial_objective = evaluate_static_objective(result.state);
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
        gradient_options.sobolev_order = sobolev_order;
        gradient_options.sobolev_cap = sobolev_cap;
        const GradientSearchResult gradient =
            active_gradient_adversary.maximize_q(
                result.state, gradient_options);
        result.state = gradient.state;
        result.initial_objective = evaluate_static_objective(result.state);
        result.evolution =
            evolve_galerkin(result.state, viscosity, final_time, dt);
        result.evaluations += gradient.trajectory_evaluations + 1;
        result.accepted_gradient_steps = gradient.accepted_steps;
    }
    result.refined_evolution =
        evolve_galerkin(result.state, viscosity, final_time, 0.5L * dt, true);
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
        StaticObjective current_objective = evaluate_static_objective(current);
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
                evaluate_static_objective(candidate);
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

using TriadKey = std::array<WaveVector, 3>;

struct InteractionAnalysis {
    std::vector<ComplexVector> advection;
    std::vector<Real> local_energy_transfer;
    std::vector<Real> nonlocal_energy_transfer;
    Real local_absolute_transfer = 0.0L;
    Real nonlocal_absolute_transfer = 0.0L;
    Real maximum_detailed_triad_residual = 0.0L;
    Real maximum_relative_detailed_triad_residual = 0.0L;
};

InteractionAnalysis analyze_interactions(const SpectralState& state) {
    InteractionAnalysis result;
    result.advection.resize(state.waves.size());
    result.local_energy_transfer.assign(state.waves.size(), 0.0L);
    result.nonlocal_energy_transfer.assign(state.waves.size(), 0.0L);
    std::map<TriadKey, std::pair<Real, Real>> detailed_triads;
    const Complex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, target_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const ComplexVector& up = state.velocity[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[target_index];
        const Complex coefficient = imaginary_unit * wave_dot(q, up);
        const ComplexVector& uq = state.velocity[q_index];
        ComplexVector& value = result.advection[target_index];
        ComplexVector pair_contribution{};
        for (std::size_t direction = 0; direction < 3; ++direction) {
            pair_contribution[direction] = coefficient * uq[direction];
            value[direction] += pair_contribution[direction];
        }

        // Since u_k is perpendicular to k, pressure projection does not
        // change its energy pairing with this interaction.
        const Real transfer = std::real(dot_hermitian(
            state.velocity[target_index], pair_contribution));
        const Integer k2 = norm_squared(k);
        const Integer p2 = norm_squared(p);
        const Integer q2 = norm_squared(q);
        const Integer smallest = std::min({k2, p2, q2});
        const Integer largest = std::max({k2, p2, q2});
        const bool local = largest <= 4 * smallest;  // max |k| / min |k| <= 2
        if (local) {
            result.local_energy_transfer[target_index] += transfer;
            result.local_absolute_transfer += std::abs(transfer);
        } else {
            result.nonlocal_energy_transfer[target_index] += transfer;
            result.nonlocal_absolute_transfer += std::abs(transfer);
        }

        TriadKey key{-k, p, q};
        std::sort(key.begin(), key.end());
        auto& [sum, absolute_sum] = detailed_triads[key];
        sum += transfer;
        absolute_sum += std::abs(transfer);
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        result.advection[index] =
            project_divergence_free(state.waves[index], result.advection[index]);
    }

    const Real total_absolute_transfer =
        result.local_absolute_transfer + result.nonlocal_absolute_transfer;
    for (const auto& [key, cancellation] : detailed_triads) {
        static_cast<void>(key);
        const auto [sum, absolute_sum] = cancellation;
        result.maximum_detailed_triad_residual =
            std::max(result.maximum_detailed_triad_residual, std::abs(sum));
        if (absolute_sum > 1e-14L * total_absolute_transfer) {
            result.maximum_relative_detailed_triad_residual =
                std::max(result.maximum_relative_detailed_triad_residual,
                         std::abs(sum) / absolute_sum);
        }
    }
    if (total_absolute_transfer > 0.0L) {
        result.maximum_detailed_triad_residual /= total_absolute_transfer;
    }
    return result;
}

struct TriadMeasurement {
    Real energy = 0.0L;
    Real enstrophy = 0.0L;
    Real palinstrophy = 0.0L;
    Real advection_norm = 0.0L;
    Real energy_pairing = 0.0L;
    Real vortex_stretching = 0.0L;
    Real divergence_residual = 0.0L;
    Real reality_residual = 0.0L;
    Real classical_ratio = 0.0L;
    Real detailed_triad_residual = 0.0L;
    Real relative_detailed_triad_residual = 0.0L;
    Real local_absolute_transfer = 0.0L;
    Real nonlocal_absolute_transfer = 0.0L;
    Real maximum_cumulative_flux = 0.0L;
    Real maximum_local_cumulative_flux = 0.0L;
    Real maximum_nonlocal_cumulative_flux = 0.0L;
    Real flux_partition_residual = 0.0L;
};

TriadMeasurement measure(const SpectralState& state) {
    const InteractionAnalysis interactions = analyze_interactions(state);
    TriadMeasurement result;
    result.detailed_triad_residual = interactions.maximum_detailed_triad_residual;
    result.relative_detailed_triad_residual =
        interactions.maximum_relative_detailed_triad_residual;
    result.local_absolute_transfer = interactions.local_absolute_transfer;
    result.nonlocal_absolute_transfer = interactions.nonlocal_absolute_transfer;
    int maximum_radius = 0;
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        const Real wave2 = static_cast<Real>(norm_squared(wave));
        maximum_radius = std::max(maximum_radius,
                                  static_cast<int>(std::ceil(std::sqrt(wave2))));
        const Real velocity2 = std::real(dot_hermitian(state.velocity[index],
                                                       state.velocity[index]));
        const Real nonlinear2 = std::real(dot_hermitian(
            interactions.advection[index], interactions.advection[index]));
        const Real pairing = std::real(dot_hermitian(
            state.velocity[index], interactions.advection[index]));
        result.energy += velocity2;
        result.enstrophy += wave2 * velocity2;
        result.palinstrophy += wave2 * wave2 * velocity2;
        result.advection_norm += nonlinear2;
        result.energy_pairing += pairing;
        result.vortex_stretching += wave2 * pairing;
        result.divergence_residual =
            std::max(result.divergence_residual, std::abs(wave_dot(wave, state.velocity[index])));
        const auto opposite = state.index.find(-wave);
        if (opposite != state.index.end()) {
            for (std::size_t direction = 0; direction < 3; ++direction) {
                result.reality_residual = std::max(
                    result.reality_residual,
                    std::abs(state.velocity[opposite->second][direction] -
                             std::conj(state.velocity[index][direction])));
            }
        }
    }
    const Real denominator = std::pow(result.enstrophy, 0.75L) *
                             std::pow(result.palinstrophy, 0.75L);
    if (denominator > 0.0L) {
        result.classical_ratio = std::abs(result.vortex_stretching) / denominator;
    }

    // Pi(K) is the nonlinear energy entering modes |k| > K. Splitting every
    // ordered triad by scale ratio gives an exact local/nonlocal partition.
    for (int cutoff = 1; cutoff < maximum_radius; ++cutoff) {
        Real local_flux = 0.0L;
        Real nonlocal_flux = 0.0L;
        const Integer cutoff2 = static_cast<Integer>(cutoff) * cutoff;
        for (std::size_t index = 0; index < state.waves.size(); ++index) {
            if (norm_squared(state.waves[index]) > cutoff2) {
                local_flux -= interactions.local_energy_transfer[index];
                nonlocal_flux -= interactions.nonlocal_energy_transfer[index];
            }
        }
        const Real total_flux = local_flux + nonlocal_flux;
        result.maximum_cumulative_flux =
            std::max(result.maximum_cumulative_flux, std::abs(total_flux));
        result.maximum_local_cumulative_flux =
            std::max(result.maximum_local_cumulative_flux, std::abs(local_flux));
        result.maximum_nonlocal_cumulative_flux =
            std::max(result.maximum_nonlocal_cumulative_flux, std::abs(nonlocal_flux));
        result.flux_partition_residual =
            std::max(result.flux_partition_residual,
                     std::abs(total_flux - local_flux - nonlocal_flux));
    }
    return result;
}

struct TriadCertificate {
    int modes = 0;
    int samples = 0;
    Real maximum_normalized_energy_residual = 0.0L;
    Real maximum_divergence_residual = 0.0L;
    Real maximum_reality_residual = 0.0L;
    Real maximum_classical_ratio = 0.0L;
    Real maximum_vortex_stretching = 0.0L;
    Real maximum_detailed_triad_residual = 0.0L;
    Real maximum_relative_detailed_triad_residual = 0.0L;
    Real maximum_nonlocal_absolute_fraction = 0.0L;
    Real maximum_flux_efficiency = 0.0L;
    Real maximum_local_cumulative_flux = 0.0L;
    Real maximum_nonlocal_cumulative_flux = 0.0L;
    Real maximum_flux_partition_residual = 0.0L;
    bool nonzero_vortex_stretching_seen = false;
};

TriadCertificate analyze_triads(int cutoff, int samples, std::uint64_t seed) {
    if (samples < 1 || samples > 100000) {
        throw std::invalid_argument("--triad-samples must be between 1 and 100000");
    }
    std::mt19937_64 generator(seed);
    TriadCertificate certificate;
    certificate.samples = samples;
    for (int sample = 0; sample < samples; ++sample) {
        const SpectralState state = SpectralStateFactory::random(cutoff, generator);
        certificate.modes = static_cast<int>(state.waves.size());
        const TriadMeasurement measurement = measure(state);
        const Real scale = std::sqrt(std::max(0.0L, measurement.energy *
                                                     measurement.advection_norm));
        const Real normalized_energy_residual =
            scale > 0.0L ? std::abs(measurement.energy_pairing) / scale : 0.0L;
        certificate.maximum_normalized_energy_residual =
            std::max(certificate.maximum_normalized_energy_residual,
                     normalized_energy_residual);
        certificate.maximum_divergence_residual =
            std::max(certificate.maximum_divergence_residual,
                     measurement.divergence_residual);
        certificate.maximum_reality_residual =
            std::max(certificate.maximum_reality_residual, measurement.reality_residual);
        certificate.maximum_classical_ratio =
            std::max(certificate.maximum_classical_ratio, measurement.classical_ratio);
        certificate.maximum_vortex_stretching =
            std::max(certificate.maximum_vortex_stretching,
                     std::abs(measurement.vortex_stretching));
        certificate.maximum_detailed_triad_residual =
            std::max(certificate.maximum_detailed_triad_residual,
                     measurement.detailed_triad_residual);
        certificate.maximum_relative_detailed_triad_residual =
            std::max(certificate.maximum_relative_detailed_triad_residual,
                     measurement.relative_detailed_triad_residual);
        const Real absolute_transfer = measurement.local_absolute_transfer +
                                       measurement.nonlocal_absolute_transfer;
        if (absolute_transfer > 0.0L) {
            certificate.maximum_nonlocal_absolute_fraction =
                std::max(certificate.maximum_nonlocal_absolute_fraction,
                         measurement.nonlocal_absolute_transfer / absolute_transfer);
            certificate.maximum_flux_efficiency =
                std::max(certificate.maximum_flux_efficiency,
                         measurement.maximum_cumulative_flux / absolute_transfer);
        }
        certificate.maximum_local_cumulative_flux =
            std::max(certificate.maximum_local_cumulative_flux,
                     measurement.maximum_local_cumulative_flux);
        certificate.maximum_nonlocal_cumulative_flux =
            std::max(certificate.maximum_nonlocal_cumulative_flux,
                     measurement.maximum_nonlocal_cumulative_flux);
        certificate.maximum_flux_partition_residual =
            std::max(certificate.maximum_flux_partition_residual,
                     measurement.flux_partition_residual);
        certificate.nonzero_vortex_stretching_seen =
            certificate.nonzero_vortex_stretching_seen ||
            std::abs(measurement.vortex_stretching) > 1e-16L;
    }
    return certificate;
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
        analyze_triads(options.triad_cutoff, options.triad_samples, options.seed);
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
    const TriadCertificate triads = analyze_triads(2, 2, 7);
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
        evolve_galerkin(adversary.state, 0.1L, 0.002L, 0.001L);
    const QDerivativeDiagnostic q_derivative =
        evaluate_q_derivative(adversary.state, 0.1L);
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
        run.objective = evaluate_static_objective(run.initial);
        run.coarse = evolve_galerkin(
            run.initial, static_cast<Real>(options.viscosity),
            static_cast<Real>(options.evolution_time),
            static_cast<Real>(options.time_step));
        run.refined = evolve_galerkin(
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
