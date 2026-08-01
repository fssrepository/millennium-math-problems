#pragma once

#include "spectral_objective.hpp"
#include "stretching_derivative_ledger.hpp"

namespace lemma {

struct LocalCriticalDerivativeSplit {
    SpectralReal nonlinear = 0.0L;
    SpectralReal viscous = 0.0L;
    SpectralReal total = 0.0L;
};

struct LocalCriticalNonlinearPartition {
    SpectralReal local = 0.0L;
    SpectralReal nonlocal = 0.0L;
    SpectralReal total = 0.0L;
};

struct LocalCriticalDerivativeLedgerReport {
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal hyperpalinstrophy = 0.0L;
    SpectralReal signed_stretching = 0.0L;
    SpectralReal global_signed_stretching = 0.0L;
    SpectralReal density = 0.0L;

    LocalCriticalDerivativeSplit stretching_derivative;
    StretchingDerivativeRoles stretching_nonlinear_roles;
    StretchingDerivativeRoles stretching_local_nonlinear_roles;
    StretchingDerivativeRoles stretching_nonlocal_nonlinear_roles;
    StretchingDerivativeRoles stretching_viscous_roles;
    LocalCriticalDerivativeSplit enstrophy_derivative;
    LocalCriticalDerivativeSplit palinstrophy_derivative;

    LocalCriticalDerivativeSplit density_from_stretching;
    LocalCriticalDerivativeSplit density_from_enstrophy;
    LocalCriticalDerivativeSplit density_from_palinstrophy;
    LocalCriticalDerivativeSplit reconstructed_density_derivative;
    LocalCriticalNonlinearPartition stretching_nonlinear_partition;
    LocalCriticalNonlinearPartition enstrophy_nonlinear_partition;
    LocalCriticalNonlinearPartition palinstrophy_nonlinear_partition;
    LocalCriticalNonlinearPartition density_nonlinear_partition;

    SpectralReal oracle_density_derivative = 0.0L;
    SpectralReal absolute_reconstruction_error = 0.0L;
    SpectralReal relative_reconstruction_error = 0.0L;
    SpectralReal enstrophy_nonlinear_identity_error = 0.0L;
    SpectralReal stretching_nonlinear_role_error = 0.0L;
    SpectralReal stretching_viscous_role_error = 0.0L;
    SpectralReal nonlinear_partition_error = 0.0L;
    SpectralReal stretching_viscous_advected_cancellation = 0.0L;
    SpectralReal enstrophy_viscous_identity_error = 0.0L;
    SpectralReal palinstrophy_viscous_identity_error = 0.0L;
    bool finite = false;
};

class LocalCriticalDerivativeLedger {
public:
    // Decomposes d[S^4/(Z P^3)]/dt into its exact S, Z, and P
    // chain-rule terms and then into Euler and viscous RHS contributions.
    [[nodiscard]] static LocalCriticalDerivativeLedgerReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralObjective& objective,
        const SpectralState& state,
        SpectralReal viscosity,
        TriadSelection selection = TriadPartition::local,
        int threads = 1);
};

}  // namespace lemma
