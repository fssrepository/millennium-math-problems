#pragma once

#include "spectral_dynamics.hpp"

namespace lemma {

struct StretchingDerivativeRoles {
    SpectralReal outer_state = 0.0L;
    SpectralReal advecting_slot = 0.0L;
    SpectralReal advected_slot = 0.0L;
    SpectralReal total = 0.0L;
};

class StretchingDerivativeLedger {
public:
    // For S=<A u,B(u,u)>, returns the three exact terms in
    // S'[h]=<A h,B(u,u)>+<A u,B(h,u)>+<A u,B(u,h)>.
    [[nodiscard]] static StretchingDerivativeRoles evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        const SpectralIncrement& direction,
        TriadSelection selection = TriadPartition::local,
        int threads = 1);
};

}  // namespace lemma
