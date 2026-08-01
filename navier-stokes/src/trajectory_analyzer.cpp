#include "trajectory_analyzer.hpp"

#include "triad_ledger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using RealVector = std::array<SpectralReal, 3>;

SpectralReal critical_integrand_from_stretching(
    SpectralReal stretching, SpectralReal enstrophy,
    SpectralReal palinstrophy) {
    const SpectralReal denominator =
        enstrophy * std::pow(palinstrophy, 3.0L);
    if (!(denominator > 0.0L)) {
        return 0.0L;
    }
    return std::pow(std::abs(stretching), 4.0L) / denominator;
}

SpectralReal energy_level_quantity_from_stretching(
    SpectralReal stretching, SpectralReal enstrophy,
    SpectralReal palinstrophy) {
    const SpectralReal denominator =
        std::pow(enstrophy, 2.0L) * std::pow(palinstrophy, 3.0L);
    if (!(denominator > 0.0L)) {
        return 0.0L;
    }
    return std::pow(std::abs(stretching), 4.0L) / denominator;
}

SpectralReal real_vector_norm(const RealVector& vector) {
    return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] +
                     vector[2] * vector[2]);
}

}  // namespace

struct TrajectoryAnalyzer::VortexPartition {
    SpectralReal local = 0.0L;
    SpectralReal nonlocal = 0.0L;
    SpectralReal near_nonlocal = 0.0L;
    SpectralReal far_nonlocal = 0.0L;
    SpectralReal selected_gap_tail = 0.0L;
    SpectralReal absolute_local_pairs = 0.0L;
    SpectralReal absolute_nonlocal_pairs = 0.0L;
};

struct TrajectoryAnalyzer::GeometryDiagnostic {
    SpectralReal maximum_vorticity = 0.0L;
    SpectralReal rms_vorticity = 0.0L;
    SpectralReal maximum_holder_half_coherence = 0.0L;
    SpectralReal mean_holder_half_coherence = 0.0L;
    SpectralReal positive_stretching_fraction = 0.0L;
    SpectralReal maximum_stretch_alignment = 0.0L;
    int high_vorticity_points = 0;
    int coherence_pairs = 0;
};

TrajectoryAnalyzer::TrajectoryAnalyzer(
    const SpectralGalerkin& configuration,
    const SpectralDynamics& dynamics,
    const SpectralObjective& objective)
    : configuration_(configuration), dynamics_(dynamics), objective_(objective) {}

StaticObjective TrajectoryAnalyzer::evaluate_static(
    const SpectralState& state) const {
    return objective_.evaluate(state);
}

TrajectoryAnalyzer::VortexPartition
TrajectoryAnalyzer::evaluate_vortex_partition(
    const SpectralState& state,
    int minimum_selected_dyadic_gap) const {
    VortexPartition result;
    const TriadLedgerReport ledger = TriadLedger::analyze(state);
    result.local = ledger.signed_local;
    result.nonlocal = ledger.signed_nonlocal;
    result.absolute_local_pairs = ledger.gaps.empty()
        ? 0.0L
        : ledger.gaps.front().absolute_pair_stretching;
    result.absolute_nonlocal_pairs =
        ledger.absolute_pair_total - result.absolute_local_pairs;
    for (const TriadGapLedgerRow& row : ledger.gaps) {
        if (row.dyadic_gap == 1) {
            result.near_nonlocal += row.signed_stretching;
        }
        if (row.dyadic_gap >= 2) {
            result.far_nonlocal += row.signed_stretching;
        }
        if (row.dyadic_gap >= minimum_selected_dyadic_gap) {
            result.selected_gap_tail += row.signed_stretching;
        }
    }
    return result;
}

TrajectoryAnalyzer::GeometryDiagnostic TrajectoryAnalyzer::evaluate_geometry(
    const SpectralState& state,
    SpectralReal high_vorticity_fraction) const {
    const int cutoff = SpectralStateOps::cutoff(state);
    const int grid_n = 2 * cutoff + 3;
    const SpectralReal spacing =
        2.0L * std::acos(-1.0L) / static_cast<SpectralReal>(grid_n);
    const std::size_t cells = static_cast<std::size_t>(grid_n) *
                              static_cast<std::size_t>(grid_n) *
                              static_cast<std::size_t>(grid_n);
    std::vector<RealVector> vorticity(cells);
    std::vector<SpectralReal> magnitude(cells, 0.0L);
    GeometryDiagnostic result;

    const auto cell_index = [grid_n](int x, int y, int z) {
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

    SpectralReal vorticity2_sum = 0.0L;
    SpectralReal positive_stretching = 0.0L;
    SpectralReal absolute_stretching = 0.0L;
    SpectralReal maximum_vorticity = 0.0L;
    SpectralReal maximum_stretch_alignment = 0.0L;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for collapse(3) num_threads(configuration_.compute_threads()) schedule(static) \
    reduction(+ : vorticity2_sum, positive_stretching, absolute_stretching) \
    reduction(max : maximum_vorticity, maximum_stretch_alignment)
#endif
    for (int z_index = 0; z_index < grid_n; ++z_index) {
        for (int y_index = 0; y_index < grid_n; ++y_index) {
            for (int x_index = 0; x_index < grid_n; ++x_index) {
                std::array<std::array<SpectralReal, 3>, 3> gradient{};
                const SpectralReal x =
                    spacing * static_cast<SpectralReal>(x_index);
                const SpectralReal y =
                    spacing * static_cast<SpectralReal>(y_index);
                const SpectralReal z =
                    spacing * static_cast<SpectralReal>(z_index);
                for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
                    const WaveVector wave = state.waves[mode];
                    const SpectralReal angle =
                        static_cast<SpectralReal>(wave.x) * x +
                        static_cast<SpectralReal>(wave.y) * y +
                        static_cast<SpectralReal>(wave.z) * z;
                    const SpectralComplex phase{
                        std::cos(angle), std::sin(angle)};
                    const std::array<int, 3> wave_component{
                        wave.x, wave.y, wave.z};
                    for (std::size_t velocity_component = 0;
                         velocity_component < 3; ++velocity_component) {
                        for (std::size_t derivative = 0; derivative < 3;
                             ++derivative) {
                            gradient[velocity_component][derivative] +=
                                std::real(
                                    imaginary_unit *
                                    static_cast<SpectralReal>(
                                        wave_component[derivative]) *
                                    state.velocity[mode][velocity_component] *
                                    phase);
                        }
                    }
                }
                const RealVector omega{
                    gradient[2][1] - gradient[1][2],
                    gradient[0][2] - gradient[2][0],
                    gradient[1][0] - gradient[0][1]};
                std::array<std::array<SpectralReal, 3>, 3> strain{};
                SpectralReal strain2 = 0.0L;
                for (std::size_t row = 0; row < 3; ++row) {
                    for (std::size_t column = 0; column < 3; ++column) {
                        strain[row][column] =
                            0.5L * (gradient[row][column] +
                                    gradient[column][row]);
                        strain2 += strain[row][column] * strain[row][column];
                    }
                }
                SpectralReal stretch = 0.0L;
                for (std::size_t row = 0; row < 3; ++row) {
                    for (std::size_t column = 0; column < 3; ++column) {
                        stretch += omega[row] * strain[row][column] *
                                   omega[column];
                    }
                }
                const auto cell = cell_index(x_index, y_index, z_index);
                vorticity[cell] = omega;
                magnitude[cell] = real_vector_norm(omega);
                maximum_vorticity =
                    std::max(maximum_vorticity, magnitude[cell]);
                vorticity2_sum += magnitude[cell] * magnitude[cell];
                positive_stretching += std::max(0.0L, stretch);
                absolute_stretching += std::abs(stretch);
                if (magnitude[cell] > 0.0L && strain2 > 0.0L) {
                    maximum_stretch_alignment = std::max(
                        maximum_stretch_alignment,
                        std::abs(stretch) /
                            (magnitude[cell] * magnitude[cell] *
                             std::sqrt(strain2)));
                }
            }
        }
    }
    result.maximum_vorticity = maximum_vorticity;
    result.maximum_stretch_alignment = maximum_stretch_alignment;
    result.rms_vorticity =
        std::sqrt(vorticity2_sum / static_cast<SpectralReal>(cells));
    if (absolute_stretching > 0.0L) {
        result.positive_stretching_fraction =
            positive_stretching / absolute_stretching;
    }

    const SpectralReal threshold =
        high_vorticity_fraction * result.maximum_vorticity;
    SpectralReal coherence_sum = 0.0L;
    SpectralReal maximum_holder_half = 0.0L;
    int high_vorticity_points = 0;
    int coherence_pairs = 0;
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for collapse(3) num_threads(configuration_.compute_threads()) schedule(static) \
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
                    const auto neighbor = cell_index(
                        neighbor_coordinates[0], neighbor_coordinates[1],
                        neighbor_coordinates[2]);
                    if (magnitude[neighbor] < threshold ||
                        magnitude[neighbor] <= 0.0L) {
                        continue;
                    }
                    const RealVector& a = vorticity[cell];
                    const RealVector& b = vorticity[neighbor];
                    const RealVector cross{
                        a[1] * b[2] - a[2] * b[1],
                        a[2] * b[0] - a[0] * b[2],
                        a[0] * b[1] - a[1] * b[0]};
                    const SpectralReal sine = std::min(
                        1.0L, real_vector_norm(cross) /
                                  (magnitude[cell] * magnitude[neighbor]));
                    const SpectralReal holder = sine / std::sqrt(spacing);
                    maximum_holder_half =
                        std::max(maximum_holder_half, holder);
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
            coherence_sum /
            static_cast<SpectralReal>(result.coherence_pairs);
    }
    return result;
}

QDerivativeDiagnostic TrajectoryAnalyzer::evaluate_q_derivative(
    const SpectralState& state, SpectralReal viscosity) const {
    const StaticObjective center = evaluate_static(state);
    const SpectralIncrement rhs = dynamics_.rhs(state, viscosity);
    SpectralReal rhs_norm2 = 0.0L;
    for (const ComplexVector& value : rhs) {
        rhs_norm2 += std::real(dot_hermitian(value, value));
    }
    QDerivativeDiagnostic result;
    if (!(center.energy > 0.0L) || !(rhs_norm2 > 0.0L) ||
        !(center.energy_level_quantity > 1e-30L)) {
        return result;
    }
    const SpectralReal base_step =
        1e-4L * std::sqrt(center.energy / rhs_norm2);
    const auto q_at_offset = [&](SpectralReal offset) {
        SpectralState shifted = state;
        for (std::size_t mode = 0; mode < shifted.velocity.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                shifted.velocity[mode][component] +=
                    offset * rhs[mode][component];
            }
        }
        return evaluate_static(shifted).energy_level_quantity;
    };
    const auto central_derivative = [&](SpectralReal step) {
        return (q_at_offset(step) - q_at_offset(-step)) / (2.0L * step);
    };
    const SpectralReal coarse = central_derivative(base_step);
    const SpectralReal refined = central_derivative(0.5L * base_step);
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

EvolutionResult TrajectoryAnalyzer::evolve(
    SpectralState state, SpectralReal viscosity,
    SpectralReal final_time, SpectralReal requested_dt,
    bool collect_vortex_partition,
    int minimum_selected_dyadic_gap) const {
    if (!(viscosity > 0.0L) || !(final_time > 0.0L) ||
        !(requested_dt > 0.0L)) {
        throw std::invalid_argument(
            "evolution viscosity, time, and dt must be positive");
    }
    EvolutionResult result;
    StaticObjective before = evaluate_static(state);
    result.initial_energy = before.energy;
    result.initial_enstrophy = before.enstrophy;
    result.initial_energy_level_quantity = before.energy_level_quantity;
    result.initial_critical_integrand = before.critical_integrand;
    result.final_critical_integrand = before.critical_integrand;
    result.maximum_energy_level_quantity = before.energy_level_quantity;
    result.maximum_critical_integrand = before.critical_integrand;
    result.maximum_enstrophy = before.enstrophy;
    const SpectralReal initial_frequency = std::sqrt(
        before.enstrophy / std::max(1e-30L, before.energy));
    const auto sample_q_derivative = [&] (
        const SpectralState& sample_state,
        const StaticObjective& sample_objective) {
        const QDerivativeDiagnostic derivative =
            evaluate_q_derivative(sample_state, viscosity);
        if (!derivative.valid) {
            return;
        }
        const SpectralReal denominator =
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
        partition_before = evaluate_vortex_partition(
            state, minimum_selected_dyadic_gap);
        result.maximum_local_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.local, before.enstrophy,
                before.palinstrophy);
        result.initial_local_critical_integrand =
            critical_integrand_from_stretching(
                partition_before.local, before.enstrophy,
                before.palinstrophy);
        result.final_local_critical_integrand =
            result.initial_local_critical_integrand;
        result.maximum_nonlocal_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.nonlocal, before.enstrophy,
                before.palinstrophy);
        result.maximum_near_nonlocal_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.near_nonlocal, before.enstrophy,
                before.palinstrophy);
        result.maximum_far_nonlocal_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.far_nonlocal, before.enstrophy,
                before.palinstrophy);
        result.maximum_selected_gap_tail_energy_level_quantity =
            energy_level_quantity_from_stretching(
                partition_before.selected_gap_tail, before.enstrophy,
                before.palinstrophy);
        result.maximum_vortex_partition_residual =
            std::abs(before.vortex_stretching -
                     std::abs(partition_before.local +
                              partition_before.nonlocal));
        const GeometryDiagnostic geometry = evaluate_geometry(state);
        result.maximum_vorticity_linf = geometry.maximum_vorticity;
        result.maximum_holder_half_coherence =
            geometry.maximum_holder_half_coherence;
        result.maximum_stretch_alignment =
            geometry.maximum_stretch_alignment;
        ++result.geometry_samples;
        sample_q_derivative(state, before);
    }

    const int cutoff = SpectralStateOps::cutoff(state);
    const SpectralReal maximum_wave2 =
        3.0L * static_cast<SpectralReal>(cutoff * cutoff);
    const SpectralReal diffusion_dt =
        0.5L / (viscosity * maximum_wave2);
    SpectralReal time = 0.0L;
    while (time < final_time) {
        const SpectralReal dt =
            std::min({requested_dt, diffusion_dt, final_time - time});
        dynamics_.rk4_step(state, viscosity, dt);
        const StaticObjective after = evaluate_static(state);
        if (collect_vortex_partition) {
            const VortexPartition partition_after =
                evaluate_vortex_partition(
                    state, minimum_selected_dyadic_gap);
            const SpectralReal local_critical_before =
                critical_integrand_from_stretching(
                    partition_before.local, before.enstrophy,
                    before.palinstrophy);
            const SpectralReal local_critical_after =
                critical_integrand_from_stretching(
                    partition_after.local, after.enstrophy,
                    after.palinstrophy);
            result.final_local_critical_integrand =
                local_critical_after;
            const SpectralReal nonlocal_critical_before =
                critical_integrand_from_stretching(
                    partition_before.nonlocal, before.enstrophy,
                    before.palinstrophy);
            const SpectralReal nonlocal_critical_after =
                critical_integrand_from_stretching(
                    partition_after.nonlocal, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal near_nonlocal_critical_before =
                critical_integrand_from_stretching(
                    partition_before.near_nonlocal, before.enstrophy,
                    before.palinstrophy);
            const SpectralReal near_nonlocal_critical_after =
                critical_integrand_from_stretching(
                    partition_after.near_nonlocal, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal far_nonlocal_critical_before =
                critical_integrand_from_stretching(
                    partition_before.far_nonlocal, before.enstrophy,
                    before.palinstrophy);
            const SpectralReal far_nonlocal_critical_after =
                critical_integrand_from_stretching(
                    partition_after.far_nonlocal, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal selected_gap_tail_critical_before =
                critical_integrand_from_stretching(
                    partition_before.selected_gap_tail, before.enstrophy,
                    before.palinstrophy);
            const SpectralReal selected_gap_tail_critical_after =
                critical_integrand_from_stretching(
                    partition_after.selected_gap_tail, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal local_q_after =
                energy_level_quantity_from_stretching(
                    partition_after.local, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal nonlocal_q_after =
                energy_level_quantity_from_stretching(
                    partition_after.nonlocal, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal near_nonlocal_q_after =
                energy_level_quantity_from_stretching(
                    partition_after.near_nonlocal, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal far_nonlocal_q_after =
                energy_level_quantity_from_stretching(
                    partition_after.far_nonlocal, after.enstrophy,
                    after.palinstrophy);
            const SpectralReal selected_gap_tail_q_after =
                energy_level_quantity_from_stretching(
                    partition_after.selected_gap_tail, after.enstrophy,
                    after.palinstrophy);
            result.integral_local_critical +=
                0.5L * dt *
                (local_critical_before + local_critical_after);
            result.integral_nonlocal_critical +=
                0.5L * dt *
                (nonlocal_critical_before + nonlocal_critical_after);
            result.integral_near_nonlocal_critical +=
                0.5L * dt *
                (near_nonlocal_critical_before +
                 near_nonlocal_critical_after);
            result.integral_far_nonlocal_critical +=
                0.5L * dt *
                (far_nonlocal_critical_before +
                 far_nonlocal_critical_after);
            result.integral_selected_gap_tail_critical +=
                0.5L * dt *
                (selected_gap_tail_critical_before +
                 selected_gap_tail_critical_after);
            result.maximum_local_energy_level_quantity = std::max(
                result.maximum_local_energy_level_quantity, local_q_after);
            result.maximum_nonlocal_energy_level_quantity = std::max(
                result.maximum_nonlocal_energy_level_quantity,
                nonlocal_q_after);
            result.maximum_near_nonlocal_energy_level_quantity = std::max(
                result.maximum_near_nonlocal_energy_level_quantity,
                near_nonlocal_q_after);
            result.maximum_far_nonlocal_energy_level_quantity = std::max(
                result.maximum_far_nonlocal_energy_level_quantity,
                far_nonlocal_q_after);
            result.maximum_selected_gap_tail_energy_level_quantity = std::max(
                result.maximum_selected_gap_tail_energy_level_quantity,
                selected_gap_tail_q_after);
            result.integral_absolute_local_vortex +=
                0.5L * dt * (std::abs(partition_before.local) +
                             std::abs(partition_after.local));
            result.integral_absolute_nonlocal_vortex +=
                0.5L * dt * (std::abs(partition_before.nonlocal) +
                             std::abs(partition_after.nonlocal));
            result.integral_absolute_total_vortex +=
                0.5L * dt *
                (std::abs(partition_before.local +
                          partition_before.nonlocal) +
                 std::abs(partition_after.local +
                          partition_after.nonlocal));
            result.maximum_vortex_partition_residual = std::max(
                result.maximum_vortex_partition_residual,
                std::abs(after.vortex_stretching -
                         std::abs(partition_after.local +
                                  partition_after.nonlocal)));
            partition_before = partition_after;
            if (result.steps % 10 == 0 || time + dt >= final_time) {
                const GeometryDiagnostic geometry = evaluate_geometry(state);
                result.maximum_vorticity_linf = std::max(
                    result.maximum_vorticity_linf,
                    geometry.maximum_vorticity);
                result.maximum_holder_half_coherence = std::max(
                    result.maximum_holder_half_coherence,
                    geometry.maximum_holder_half_coherence);
                result.maximum_stretch_alignment = std::max(
                    result.maximum_stretch_alignment,
                    geometry.maximum_stretch_alignment);
                ++result.geometry_samples;
                sample_q_derivative(state, after);
            }
        }
        result.integral_critical +=
            0.5L * dt *
            (before.critical_integrand + after.critical_integrand);
        result.integral_enstrophy +=
            0.5L * dt * (before.enstrophy + after.enstrophy);
        result.maximum_energy_level_quantity = std::max(
            result.maximum_energy_level_quantity,
            after.energy_level_quantity);
        result.maximum_critical_integrand = std::max(
            result.maximum_critical_integrand, after.critical_integrand);
        result.final_critical_integrand = after.critical_integrand;
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
    result.energy_balance_residual =
        result.final_energy - result.initial_energy +
        2.0L * viscosity * result.integral_enstrophy;
    return result;
}

}  // namespace lemma
