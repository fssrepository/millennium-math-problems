#pragma once

#include "spectral_objective.hpp"

namespace lemma {

struct QDerivativeDiagnostic {
    SpectralReal derivative = 0.0L;
    SpectralReal log_derivative = 0.0L;
    SpectralReal relative_refinement_error = 0.0L;
    bool valid = false;
};

struct EvolutionResult {
    int steps = 0;
    SpectralReal time = 0.0L;
    SpectralReal initial_energy = 0.0L;
    SpectralReal final_energy = 0.0L;
    SpectralReal initial_enstrophy = 0.0L;
    SpectralReal final_enstrophy = 0.0L;
    SpectralReal initial_energy_level_quantity = 0.0L;
    SpectralReal final_energy_level_quantity = 0.0L;
    SpectralReal initial_critical_integrand = 0.0L;
    SpectralReal final_critical_integrand = 0.0L;
    SpectralReal initial_local_critical_integrand = 0.0L;
    SpectralReal final_local_critical_integrand = 0.0L;
    SpectralReal integral_critical = 0.0L;
    SpectralReal integral_enstrophy = 0.0L;
    SpectralReal maximum_energy_level_quantity = 0.0L;
    SpectralReal maximum_critical_integrand = 0.0L;
    SpectralReal maximum_enstrophy = 0.0L;
    SpectralReal energy_balance_residual = 0.0L;
    SpectralReal integral_absolute_local_vortex = 0.0L;
    SpectralReal integral_absolute_nonlocal_vortex = 0.0L;
    SpectralReal integral_absolute_total_vortex = 0.0L;
    SpectralReal integral_local_critical = 0.0L;
    SpectralReal integral_nonlocal_critical = 0.0L;
    SpectralReal integral_near_nonlocal_critical = 0.0L;
    SpectralReal integral_far_nonlocal_critical = 0.0L;
    SpectralReal integral_selected_gap_tail_critical = 0.0L;
    SpectralReal maximum_local_energy_level_quantity = 0.0L;
    SpectralReal maximum_nonlocal_energy_level_quantity = 0.0L;
    SpectralReal maximum_near_nonlocal_energy_level_quantity = 0.0L;
    SpectralReal maximum_far_nonlocal_energy_level_quantity = 0.0L;
    SpectralReal maximum_selected_gap_tail_energy_level_quantity = 0.0L;
    SpectralReal maximum_vortex_partition_residual = 0.0L;
    SpectralReal maximum_vorticity_linf = 0.0L;
    SpectralReal maximum_holder_half_coherence = 0.0L;
    SpectralReal maximum_stretch_alignment = 0.0L;
    SpectralReal maximum_positive_q_log_growth_ratio = 0.0L;
    SpectralReal maximum_q_derivative_refinement_error = 0.0L;
    int geometry_samples = 0;
    int q_derivative_samples = 0;
    bool finite = true;
};

class TrajectoryAnalyzer {
public:
    TrajectoryAnalyzer(const SpectralGalerkin& configuration,
                       const SpectralDynamics& dynamics,
                       const SpectralObjective& objective);

    [[nodiscard]] StaticObjective evaluate_static(
        const SpectralState& state) const;
    [[nodiscard]] QDerivativeDiagnostic evaluate_q_derivative(
        const SpectralState& state, SpectralReal viscosity) const;
    [[nodiscard]] EvolutionResult evolve(
        SpectralState state, SpectralReal viscosity,
        SpectralReal final_time, SpectralReal requested_dt,
        bool collect_vortex_partition = false,
        int minimum_selected_dyadic_gap = 2) const;

private:
    struct VortexPartition;
    struct GeometryDiagnostic;

    [[nodiscard]] VortexPartition evaluate_vortex_partition(
        const SpectralState& state,
        int minimum_selected_dyadic_gap) const;
    [[nodiscard]] GeometryDiagnostic evaluate_geometry(
        const SpectralState& state,
        SpectralReal high_vorticity_fraction = 0.5L) const;

    const SpectralGalerkin& configuration_;
    const SpectralDynamics& dynamics_;
    const SpectralObjective& objective_;
};

}  // namespace lemma
