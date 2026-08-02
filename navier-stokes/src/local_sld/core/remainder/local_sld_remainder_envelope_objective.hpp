#pragma once

#include "spectral_dynamics.hpp"

namespace lemma {

struct LocalSldRemainderEnvelopeValue {
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal stretching = 0.0L;
    SpectralReal projection_coefficient = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal projected_h3_correction = 0.0L;
    SpectralReal projected_commutator_pairing = 0.0L;
    SpectralReal commutator_hminus1_norm2 = 0.0L;
    SpectralReal upper_envelope = 0.0L;
    SpectralReal target_scale = 0.0L;
    SpectralReal target_ratio = 0.0L;
    bool finite = false;
};

class LocalSldRemainderEnvelopeObjective {
public:
    explicit LocalSldRemainderEnvelopeObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection =
            TriadSelection::local_without_equal_low_doubling());

    [[nodiscard]] LocalSldRemainderEnvelopeValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement target_ratio_gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
};

}  // namespace lemma
