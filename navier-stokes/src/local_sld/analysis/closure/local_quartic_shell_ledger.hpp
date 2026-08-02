#pragma once

#include "local_critical_derivative_ledger.hpp"
#include "shifted_critical_density.hpp"

#include <vector>

namespace lemma {

struct LocalQuarticShellRow {
    int shell = 0;
    int modes = 0;
    SpectralReal outer_state = 0.0L;
    SpectralReal advecting_slot = 0.0L;
    SpectralReal advected_slot = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal normalized_outer_state = 0.0L;
    SpectralReal normalized_advecting_slot = 0.0L;
    SpectralReal normalized_advected_slot = 0.0L;
    SpectralReal normalized_enstrophy = 0.0L;
    SpectralReal normalized_palinstrophy = 0.0L;
    SpectralReal normalized_total = 0.0L;
};

struct LocalQuarticEigenShellRow {
    SpectralInteger wave_squared = 0;
    int modes = 0;
    SpectralReal normalized_outer_state = 0.0L;
    SpectralReal normalized_advecting_slot = 0.0L;
    SpectralReal normalized_advected_slot = 0.0L;
    SpectralReal normalized_enstrophy = 0.0L;
    SpectralReal normalized_palinstrophy = 0.0L;
    SpectralReal normalized_total = 0.0L;
};

struct LocalQuarticDyadicShellRow {
    int shell = 0;
    int modes = 0;
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal local_advection_h1_squared = 0.0L;
    SpectralReal normalized_outer_state = 0.0L;
    SpectralReal normalized_advecting_slot = 0.0L;
    SpectralReal normalized_advected_slot = 0.0L;
    SpectralReal normalized_enstrophy = 0.0L;
    SpectralReal normalized_palinstrophy = 0.0L;
    SpectralReal normalized_total = 0.0L;
};

struct LocalQuarticShellReport {
    std::vector<LocalQuarticShellRow> shells;
    std::vector<LocalQuarticEigenShellRow> eigen_shells;
    std::vector<LocalQuarticDyadicShellRow> dyadic_shells;
    SpectralReal normalized_total = 0.0L;
    SpectralReal positive_dyadic_shell_total = 0.0L;
    SpectralReal negative_dyadic_shell_total = 0.0L;
    SpectralReal within_target_mode_cancellation_fraction = 0.0L;
    SpectralReal between_modes_within_shell_cancellation_fraction = 0.0L;
    SpectralReal between_modes_within_eigen_shell_cancellation_fraction =
        0.0L;
    SpectralReal between_eigen_shells_cancellation_fraction = 0.0L;
    SpectralReal between_modes_within_dyadic_shell_cancellation_fraction =
        0.0L;
    SpectralReal between_dyadic_shells_cancellation_fraction = 0.0L;
    SpectralReal within_shell_cancellation_fraction = 0.0L;
    SpectralReal between_shell_cancellation_fraction = 0.0L;
    SpectralReal raw_role_reconstruction_error = 0.0L;
    SpectralReal normalized_reconstruction_error = 0.0L;
    bool finite = false;
};

class LocalQuarticShellLedger {
public:
    // Groups the complete local-local quartic derivative by target hard shell.
    // This distinguishes cancellation inside one output scale from cancellation
    // that occurs only after unrelated shells are summed.
    [[nodiscard]] static LocalQuarticShellReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        const ShiftedCriticalDensityDiagnostic& diagnostic,
        const LocalCriticalDerivativeLedgerReport& derivative_ledger);
};

}  // namespace lemma
