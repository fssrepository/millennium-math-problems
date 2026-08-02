#include "local_quartic_shell_ledger.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace lemma {
namespace {

SpectralReal relative_error(SpectralReal value, SpectralReal reference) {
    const SpectralReal scale = std::max(
        {std::abs(value), std::abs(reference), 1.0e-30L});
    return std::abs(value - reference) / scale;
}

int hard_shell(WaveVector wave) {
    return std::max(
        {std::abs(wave.x), std::abs(wave.y), std::abs(wave.z)});
}

int dyadic_shell(SpectralInteger wave_squared) {
    int shell = 0;
    SpectralInteger next_shell_squared = 4;
    while (wave_squared >= next_shell_squared) {
        ++shell;
        next_shell_squared *= 4;
    }
    return shell;
}

}  // namespace

LocalQuarticShellReport LocalQuarticShellLedger::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const ShiftedCriticalDensityDiagnostic& diagnostic,
    const LocalCriticalDerivativeLedgerReport& derivative_ledger) {
    const int cutoff = SpectralStateOps::cutoff(state);
    LocalQuarticShellReport result;
    result.shells.resize(static_cast<std::size_t>(cutoff + 1));
    for (int shell = 0; shell <= cutoff; ++shell) {
        result.shells[static_cast<std::size_t>(shell)].shell = shell;
    }
    const SpectralInteger maximum_wave_squared =
        static_cast<SpectralInteger>(3) * cutoff * cutoff;
    result.eigen_shells.resize(
        static_cast<std::size_t>(maximum_wave_squared + 1));
    for (SpectralInteger wave_squared = 0;
         wave_squared <= maximum_wave_squared; ++wave_squared) {
        result.eigen_shells[static_cast<std::size_t>(wave_squared)]
            .wave_squared = wave_squared;
    }
    int maximum_dyadic_shell = 0;
    while ((static_cast<SpectralInteger>(1) <<
            (2 * (maximum_dyadic_shell + 1))) <=
           maximum_wave_squared) {
        ++maximum_dyadic_shell;
    }
    result.dyadic_shells.resize(
        static_cast<std::size_t>(maximum_dyadic_shell + 1));
    for (int shell = 0; shell <= maximum_dyadic_shell; ++shell) {
        result.dyadic_shells[static_cast<std::size_t>(shell)].shell = shell;
    }

    const SpectralIncrement local_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::local);
    SpectralIncrement local_direction = local_advection;
    for (ComplexVector& mode : local_direction) {
        for (SpectralComplex& component : mode) {
            component = -component;
        }
    }
    const SpectralIncrement advecting_variation =
        dynamics.advection_bilinear_direct_partition(
            state, local_direction, state.velocity,
            TriadPartition::local);
    const SpectralIncrement advected_variation =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, local_direction,
            TriadPartition::local);

    const SpectralReal shifted_density =
        diagnostic.local_critical_density + diagnostic.initial_ep_shift;
    const SpectralReal normalization_scale =
        shifted_density * diagnostic.normalization;
    SpectralReal stretching_coefficient = 0.0L;
    SpectralReal enstrophy_coefficient = 0.0L;
    SpectralReal palinstrophy_coefficient = 0.0L;
    if (normalization_scale > 0.0L &&
        derivative_ledger.enstrophy > 0.0L &&
        derivative_ledger.palinstrophy > 0.0L) {
        const SpectralReal stretching2 =
            derivative_ledger.signed_stretching *
            derivative_ledger.signed_stretching;
        const SpectralReal stretching3 = stretching2 *
            derivative_ledger.signed_stretching;
        const SpectralReal stretching4 = stretching2 * stretching2;
        const SpectralReal palinstrophy3 =
            derivative_ledger.palinstrophy *
            derivative_ledger.palinstrophy *
            derivative_ledger.palinstrophy;
        stretching_coefficient = 4.0L * stretching3 /
            (derivative_ledger.enstrophy * palinstrophy3 *
             normalization_scale);
        enstrophy_coefficient = -stretching4 /
            (derivative_ledger.enstrophy *
             derivative_ledger.enstrophy * palinstrophy3 *
             normalization_scale);
        palinstrophy_coefficient = -3.0L * stretching4 /
            (derivative_ledger.enstrophy * palinstrophy3 *
             derivative_ledger.palinstrophy *
             normalization_scale);
    }

    SpectralReal mode_component_absolute_sum = 0.0L;
    SpectralReal mode_total_absolute_sum = 0.0L;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        LocalQuarticShellRow& row = result.shells[
            static_cast<std::size_t>(hard_shell(state.waves[mode]))];
        ++row.modes;
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        const ComplexVector& velocity = state.velocity[mode];
        const ComplexVector& direction = local_direction[mode];
        const SpectralReal outer = wave2 * std::real(dot_hermitian(
            local_advection[mode], direction));
        const SpectralReal advecting = wave2 * std::real(dot_hermitian(
            velocity, advecting_variation[mode]));
        const SpectralReal advected = wave2 * std::real(dot_hermitian(
            velocity, advected_variation[mode]));
        const SpectralReal enstrophy = 2.0L * wave2 * std::real(
            dot_hermitian(velocity, direction));
        const SpectralReal palinstrophy = 2.0L * wave4 * std::real(
            dot_hermitian(velocity, direction));
        row.outer_state += outer;
        row.advecting_slot += advecting;
        row.advected_slot += advected;
        row.enstrophy += enstrophy;
        row.palinstrophy += palinstrophy;
        const SpectralReal normalized_outer =
            stretching_coefficient * outer;
        const SpectralReal normalized_advecting =
            stretching_coefficient * advecting;
        const SpectralReal normalized_advected =
            stretching_coefficient * advected;
        const SpectralReal normalized_enstrophy =
            enstrophy_coefficient * enstrophy;
        const SpectralReal normalized_palinstrophy =
            palinstrophy_coefficient * palinstrophy;
        LocalQuarticEigenShellRow& eigen_row = result.eigen_shells[
            static_cast<std::size_t>(
                norm_squared(state.waves[mode]))];
        ++eigen_row.modes;
        eigen_row.normalized_outer_state += normalized_outer;
        eigen_row.normalized_advecting_slot += normalized_advecting;
        eigen_row.normalized_advected_slot += normalized_advected;
        eigen_row.normalized_enstrophy += normalized_enstrophy;
        eigen_row.normalized_palinstrophy += normalized_palinstrophy;
        eigen_row.normalized_total +=
            normalized_outer + normalized_advecting +
            normalized_advected + normalized_enstrophy +
            normalized_palinstrophy;
        LocalQuarticDyadicShellRow& dyadic_row = result.dyadic_shells[
            static_cast<std::size_t>(dyadic_shell(
                norm_squared(state.waves[mode])))];
        ++dyadic_row.modes;
        const SpectralReal mode_energy = std::real(dot_hermitian(
            velocity, velocity));
        dyadic_row.energy += mode_energy;
        dyadic_row.enstrophy += wave2 * mode_energy;
        dyadic_row.palinstrophy += wave4 * mode_energy;
        dyadic_row.local_advection_h1_squared += -outer;
        dyadic_row.normalized_outer_state += normalized_outer;
        dyadic_row.normalized_advecting_slot += normalized_advecting;
        dyadic_row.normalized_advected_slot += normalized_advected;
        dyadic_row.normalized_enstrophy += normalized_enstrophy;
        dyadic_row.normalized_palinstrophy += normalized_palinstrophy;
        dyadic_row.normalized_total +=
            normalized_outer + normalized_advecting +
            normalized_advected + normalized_enstrophy +
            normalized_palinstrophy;
        mode_component_absolute_sum +=
            std::abs(normalized_outer) +
            std::abs(normalized_advecting) +
            std::abs(normalized_advected) +
            std::abs(normalized_enstrophy) +
            std::abs(normalized_palinstrophy);
        mode_total_absolute_sum += std::abs(
            normalized_outer + normalized_advecting +
            normalized_advected + normalized_enstrophy +
            normalized_palinstrophy);
    }

    SpectralReal raw_outer = 0.0L;
    SpectralReal raw_advecting = 0.0L;
    SpectralReal raw_advected = 0.0L;
    SpectralReal raw_enstrophy = 0.0L;
    SpectralReal raw_palinstrophy = 0.0L;
    SpectralReal component_absolute_sum = 0.0L;
    SpectralReal shell_total_absolute_sum = 0.0L;
    for (LocalQuarticShellRow& row : result.shells) {
        row.normalized_outer_state =
            stretching_coefficient * row.outer_state;
        row.normalized_advecting_slot =
            stretching_coefficient * row.advecting_slot;
        row.normalized_advected_slot =
            stretching_coefficient * row.advected_slot;
        row.normalized_enstrophy =
            enstrophy_coefficient * row.enstrophy;
        row.normalized_palinstrophy =
            palinstrophy_coefficient * row.palinstrophy;
        row.normalized_total = row.normalized_outer_state +
            row.normalized_advecting_slot +
            row.normalized_advected_slot + row.normalized_enstrophy +
            row.normalized_palinstrophy;
        raw_outer += row.outer_state;
        raw_advecting += row.advecting_slot;
        raw_advected += row.advected_slot;
        raw_enstrophy += row.enstrophy;
        raw_palinstrophy += row.palinstrophy;
        result.normalized_total += row.normalized_total;
        component_absolute_sum +=
            std::abs(row.normalized_outer_state) +
            std::abs(row.normalized_advecting_slot) +
            std::abs(row.normalized_advected_slot) +
            std::abs(row.normalized_enstrophy) +
            std::abs(row.normalized_palinstrophy);
        shell_total_absolute_sum += std::abs(row.normalized_total);
    }
    SpectralReal eigen_shell_total_absolute_sum = 0.0L;
    for (const LocalQuarticEigenShellRow& row : result.eigen_shells) {
        eigen_shell_total_absolute_sum +=
            std::abs(row.normalized_total);
    }
    SpectralReal dyadic_shell_total_absolute_sum = 0.0L;
    for (const LocalQuarticDyadicShellRow& row : result.dyadic_shells) {
        dyadic_shell_total_absolute_sum +=
            std::abs(row.normalized_total);
        result.positive_dyadic_shell_total +=
            std::max(0.0L, row.normalized_total);
        result.negative_dyadic_shell_total +=
            std::min(0.0L, row.normalized_total);
    }
    if (component_absolute_sum > 0.0L) {
        result.within_shell_cancellation_fraction =
            1.0L - shell_total_absolute_sum / component_absolute_sum;
    }
    if (mode_component_absolute_sum > 0.0L) {
        result.within_target_mode_cancellation_fraction =
            1.0L - mode_total_absolute_sum /
                mode_component_absolute_sum;
    }
    if (mode_total_absolute_sum > 0.0L) {
        result.between_modes_within_shell_cancellation_fraction =
            1.0L - shell_total_absolute_sum /
                mode_total_absolute_sum;
        result.between_modes_within_eigen_shell_cancellation_fraction =
            1.0L - eigen_shell_total_absolute_sum /
                mode_total_absolute_sum;
        result.between_modes_within_dyadic_shell_cancellation_fraction =
            1.0L - dyadic_shell_total_absolute_sum /
                mode_total_absolute_sum;
    }
    if (eigen_shell_total_absolute_sum > 0.0L) {
        result.between_eigen_shells_cancellation_fraction =
            1.0L - std::abs(result.normalized_total) /
                eigen_shell_total_absolute_sum;
    }
    if (dyadic_shell_total_absolute_sum > 0.0L) {
        result.between_dyadic_shells_cancellation_fraction =
            1.0L - std::abs(result.normalized_total) /
                dyadic_shell_total_absolute_sum;
    }
    if (shell_total_absolute_sum > 0.0L) {
        result.between_shell_cancellation_fraction =
            1.0L - std::abs(result.normalized_total) /
                shell_total_absolute_sum;
    }
    const StretchingDerivativeRoles& roles = derivative_ledger
        .stretching_local_nonlinear_roles;
    result.raw_role_reconstruction_error = std::max({
        relative_error(raw_outer, roles.outer_state),
        relative_error(raw_advecting, roles.advecting_slot),
        relative_error(raw_advected, roles.advected_slot),
        relative_error(
            raw_enstrophy,
            derivative_ledger.enstrophy_nonlinear_partition.local),
        relative_error(
            raw_palinstrophy,
            derivative_ledger.palinstrophy_nonlinear_partition.local),
    });
    const SpectralReal expected_normalized =
        derivative_ledger.density_nonlinear_partition.local /
        normalization_scale;
    result.normalized_reconstruction_error = relative_error(
        result.normalized_total, expected_normalized);
    result.finite =
        std::isfinite(result.normalized_total) &&
        result.raw_role_reconstruction_error < 1.0e-12L &&
        result.normalized_reconstruction_error < 1.0e-12L;
    return result;
}

}  // namespace lemma
