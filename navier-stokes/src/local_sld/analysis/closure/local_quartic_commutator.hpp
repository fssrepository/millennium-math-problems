#pragma once

#include "local_critical_derivative_ledger.hpp"

namespace lemma {

struct LocalQuarticCommutatorReport {
    SpectralReal outer_state_derivative = 0.0L;
    SpectralReal advected_slot_derivative = 0.0L;
    SpectralReal combined_derivative = 0.0L;
    SpectralReal negative_commutator_pairing = 0.0L;
    SpectralReal local_advection_l2_squared = 0.0L;
    SpectralReal commutator_l2_squared = 0.0L;
    SpectralReal cauchy_ratio = 0.0L;
    SpectralReal identity_relative_error = 0.0L;
    SpectralReal maximum_symbol_ratio = 0.0L;
    bool symbol_bound_holds = false;
    bool finite = false;
};

class LocalQuarticCommutator {
public:
    // Certifies
    // -<A B_L,B_L>-<A u,B_L(u,B_L)>
    // =-<B_L,A B_L-B_L(u,A u)>.
    [[nodiscard]] static LocalQuarticCommutatorReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        const LocalCriticalDerivativeLedgerReport& derivative_ledger);
};

}  // namespace lemma
