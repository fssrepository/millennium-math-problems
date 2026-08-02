#include "local_critical_derivative_ledger.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal pairing(const SpectralIncrement& gradient,
                     const SpectralIncrement& direction) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        result += std::real(
            dot_hermitian(gradient[mode], direction[mode]));
    }
    return result;
}

LocalCriticalDerivativeSplit split_pairing(
    const SpectralIncrement& gradient,
    const SpectralIncrement& nonlinear,
    const SpectralIncrement& viscous) {
    LocalCriticalDerivativeSplit result;
    result.nonlinear = pairing(gradient, nonlinear);
    result.viscous = pairing(gradient, viscous);
    result.total = result.nonlinear + result.viscous;
    return result;
}

LocalCriticalDerivativeSplit scaled(
    const LocalCriticalDerivativeSplit& source,
    SpectralReal coefficient) {
    return {
        coefficient * source.nonlinear,
        coefficient * source.viscous,
        coefficient * source.total,
    };
}

LocalCriticalDerivativeSplit add(
    const LocalCriticalDerivativeSplit& left,
    const LocalCriticalDerivativeSplit& right) {
    return {
        left.nonlinear + right.nonlinear,
        left.viscous + right.viscous,
        left.total + right.total,
    };
}

LocalCriticalNonlinearPartition partition_pairing(
    const SpectralIncrement& gradient,
    const SpectralIncrement& local,
    const SpectralIncrement& nonlocal) {
    LocalCriticalNonlinearPartition result;
    result.local = pairing(gradient, local);
    result.nonlocal = pairing(gradient, nonlocal);
    result.total = result.local + result.nonlocal;
    return result;
}

LocalCriticalNonlinearPartition scaled_partition(
    const LocalCriticalNonlinearPartition& source,
    SpectralReal coefficient) {
    return {
        coefficient * source.local,
        coefficient * source.nonlocal,
        coefficient * source.total,
    };
}

LocalCriticalNonlinearPartition add_partition(
    const LocalCriticalNonlinearPartition& left,
    const LocalCriticalNonlinearPartition& right) {
    return {
        left.local + right.local,
        left.nonlocal + right.nonlocal,
        left.total + right.total,
    };
}

SpectralReal relative_error(SpectralReal value, SpectralReal reference) {
    const SpectralReal scale = std::max(
        {std::abs(value), std::abs(reference), 1.0e-30L});
    return std::abs(value - reference) / scale;
}

}  // namespace

LocalCriticalDerivativeLedgerReport
LocalCriticalDerivativeLedger::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralObjective& objective,
    const SpectralState& state,
    SpectralReal viscosity,
    TriadSelection selection,
    int threads) {
    if (!(viscosity > 0.0L) || !std::isfinite(viscosity)) {
        throw std::invalid_argument(
            "local critical derivative viscosity must be positive");
    }

    const StaticObjective value = objective.evaluate(state, selection);
    const SpectralIncrement advection = dynamics.advection(state);
    const SpectralIncrement local_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::local);
    const SpectralIncrement nonlocal_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::nonlocal);
    SpectralIncrement nonlinear(state.waves.size());
    SpectralIncrement local_nonlinear(state.waves.size());
    SpectralIncrement nonlocal_nonlinear(state.waves.size());
    SpectralIncrement viscous(state.waves.size());
    SpectralIncrement enstrophy_gradient(state.waves.size());
    SpectralIncrement palinstrophy_gradient(state.waves.size());
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        for (std::size_t component = 0; component < 3; ++component) {
            nonlinear[mode][component] = -advection[mode][component];
            local_nonlinear[mode][component] =
                -local_advection[mode][component];
            nonlocal_nonlinear[mode][component] =
                -nonlocal_advection[mode][component];
            viscous[mode][component] =
                -viscosity * wave2 * state.velocity[mode][component];
            enstrophy_gradient[mode][component] =
                2.0L * wave2 * state.velocity[mode][component];
            palinstrophy_gradient[mode][component] =
                2.0L * wave4 * state.velocity[mode][component];
        }
    }

    const SpectralIncrement stretching_gradient =
        objective.signed_stretching_gradient(state, selection);
    const SpectralIncrement density_gradient =
        objective.critical_integrand_gradient(state, selection);

    LocalCriticalDerivativeLedgerReport result;
    result.energy = value.energy;
    result.enstrophy = value.enstrophy;
    result.palinstrophy = value.palinstrophy;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal mode_energy = std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
        result.hyperpalinstrophy +=
            wave2 * wave2 * wave2 * mode_energy;
        result.global_signed_stretching += wave2 * std::real(
            dot_hermitian(state.velocity[mode], advection[mode]));
    }
    result.signed_stretching = value.signed_vortex_stretching;
    result.density = value.critical_integrand;
    result.stretching_derivative = split_pairing(
        stretching_gradient, nonlinear, viscous);
    result.stretching_nonlinear_roles =
        StretchingDerivativeLedger::evaluate(
            dynamics, state, nonlinear, selection, threads);
    result.stretching_local_nonlinear_roles =
        StretchingDerivativeLedger::evaluate(
            dynamics, state, local_nonlinear, selection, threads);
    result.stretching_nonlocal_nonlinear_roles =
        StretchingDerivativeLedger::evaluate(
            dynamics, state, nonlocal_nonlinear, selection, threads);
    result.stretching_viscous_roles =
        StretchingDerivativeLedger::evaluate(
            dynamics, state, viscous, selection, threads);
    result.enstrophy_derivative = split_pairing(
        enstrophy_gradient, nonlinear, viscous);
    result.palinstrophy_derivative = split_pairing(
        palinstrophy_gradient, nonlinear, viscous);
    result.stretching_nonlinear_partition = partition_pairing(
        stretching_gradient, local_nonlinear, nonlocal_nonlinear);
    result.enstrophy_nonlinear_partition = partition_pairing(
        enstrophy_gradient, local_nonlinear, nonlocal_nonlinear);
    result.palinstrophy_nonlinear_partition = partition_pairing(
        palinstrophy_gradient, local_nonlinear, nonlocal_nonlinear);

    if (value.enstrophy > 0.0L && value.palinstrophy > 0.0L) {
        const SpectralReal stretching2 =
            value.signed_vortex_stretching *
            value.signed_vortex_stretching;
        const SpectralReal stretching3 =
            stretching2 * value.signed_vortex_stretching;
        const SpectralReal stretching4 = stretching2 * stretching2;
        const SpectralReal palinstrophy2 =
            value.palinstrophy * value.palinstrophy;
        const SpectralReal palinstrophy3 =
            palinstrophy2 * value.palinstrophy;
        result.density_from_stretching = scaled(
            result.stretching_derivative,
            4.0L * stretching3 /
                (value.enstrophy * palinstrophy3));
        result.density_from_enstrophy = scaled(
            result.enstrophy_derivative,
            -stretching4 /
                (value.enstrophy * value.enstrophy * palinstrophy3));
        result.density_from_palinstrophy = scaled(
            result.palinstrophy_derivative,
            -3.0L * stretching4 /
                (value.enstrophy * palinstrophy3 *
                 value.palinstrophy));
        const SpectralReal stretching_coefficient =
            4.0L * stretching3 /
            (value.enstrophy * palinstrophy3);
        const SpectralReal enstrophy_coefficient =
            -stretching4 /
            (value.enstrophy * value.enstrophy * palinstrophy3);
        const SpectralReal palinstrophy_coefficient =
            -3.0L * stretching4 /
            (value.enstrophy * palinstrophy3 * value.palinstrophy);
        result.density_nonlinear_partition = add_partition(
            add_partition(
                scaled_partition(
                    result.stretching_nonlinear_partition,
                    stretching_coefficient),
                scaled_partition(
                    result.enstrophy_nonlinear_partition,
                    enstrophy_coefficient)),
            scaled_partition(
                result.palinstrophy_nonlinear_partition,
                palinstrophy_coefficient));
    }
    result.reconstructed_density_derivative = add(
        add(result.density_from_stretching,
            result.density_from_enstrophy),
        result.density_from_palinstrophy);
    result.oracle_density_derivative =
        pairing(density_gradient, nonlinear) +
        pairing(density_gradient, viscous);
    result.absolute_reconstruction_error = std::abs(
        result.reconstructed_density_derivative.total -
        result.oracle_density_derivative);
    result.relative_reconstruction_error = relative_error(
        result.reconstructed_density_derivative.total,
        result.oracle_density_derivative);
    result.enstrophy_nonlinear_identity_error = relative_error(
        result.enstrophy_derivative.nonlinear,
        -2.0L * result.global_signed_stretching);
    result.stretching_nonlinear_role_error = relative_error(
        result.stretching_nonlinear_roles.total,
        result.stretching_derivative.nonlinear);
    result.stretching_viscous_role_error = relative_error(
        result.stretching_viscous_roles.total,
        result.stretching_derivative.viscous);
    result.nonlinear_partition_error = relative_error(
        result.density_nonlinear_partition.total,
        result.reconstructed_density_derivative.nonlinear);
    result.stretching_viscous_advected_cancellation =
        std::abs(result.stretching_viscous_roles.advected_slot) /
        std::max(
            std::abs(result.stretching_viscous_roles.outer_state) +
                std::abs(
                    result.stretching_viscous_roles.advecting_slot),
            1.0e-30L);
    result.enstrophy_viscous_identity_error = relative_error(
        result.enstrophy_derivative.viscous,
        -2.0L * viscosity * result.palinstrophy);
    result.palinstrophy_viscous_identity_error = relative_error(
        result.palinstrophy_derivative.viscous,
        -2.0L * viscosity * result.hyperpalinstrophy);
    result.finite =
        std::isfinite(result.oracle_density_derivative) &&
        std::isfinite(result.relative_reconstruction_error) &&
        result.relative_reconstruction_error < 1.0e-12L &&
        result.enstrophy_nonlinear_identity_error < 1.0e-15L &&
        result.stretching_nonlinear_role_error < 1.0e-12L &&
        result.stretching_viscous_role_error < 1.0e-12L &&
        result.nonlinear_partition_error < 1.0e-12L &&
        result.stretching_viscous_advected_cancellation < 1.0e-12L &&
        result.enstrophy_viscous_identity_error < 1.0e-15L &&
        result.palinstrophy_viscous_identity_error < 1.0e-15L;
    return result;
}

}  // namespace lemma
