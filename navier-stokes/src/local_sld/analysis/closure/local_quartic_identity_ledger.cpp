#include "local_quartic_identity_ledger.hpp"

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

}  // namespace

LocalQuarticIdentityReport LocalQuarticIdentityLedger::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const LocalCriticalDerivativeLedgerReport& derivative_ledger) {
    const SpectralIncrement local_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::local);
    LocalQuarticIdentityReport result;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        result.local_advection_h1_squared += wave2 * std::real(
            dot_hermitian(
                local_advection[mode], local_advection[mode]));
    }
    result.local_outer_state_derivative = derivative_ledger
        .stretching_local_nonlinear_roles.outer_state;
    result.local_outer_negative_square_error = relative_error(
        result.local_outer_state_derivative,
        -result.local_advection_h1_squared);

    result.local_enstrophy_derivative = derivative_ledger
        .enstrophy_nonlinear_partition.local;
    result.local_enstrophy_identity_error = relative_error(
        result.local_enstrophy_derivative,
        -2.0L * derivative_ledger.signed_stretching);

    result.nonlocal_enstrophy_derivative = derivative_ledger
        .enstrophy_nonlinear_partition.nonlocal;
    const SpectralReal nonlocal_signed_stretching =
        derivative_ledger.global_signed_stretching -
        derivative_ledger.signed_stretching;
    result.nonlocal_enstrophy_identity_error = relative_error(
        result.nonlocal_enstrophy_derivative,
        -2.0L * nonlocal_signed_stretching);
    result.finite =
        std::isfinite(result.local_advection_h1_squared) &&
        result.local_outer_negative_square_error < 1.0e-15L &&
        result.local_enstrophy_identity_error < 1.0e-15L &&
        result.nonlocal_enstrophy_identity_error < 1.0e-15L;
    return result;
}

}  // namespace lemma
