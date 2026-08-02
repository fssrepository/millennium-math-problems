#pragma once

#include "local_quartic_closure_objective.hpp"

namespace lemma {

struct LocalSldRemainderAbsorptionValue {
    SpectralReal theta = 1.0L;
    SpectralReal signed_bracket = 0.0L;
    SpectralReal first_square_norm2 = 0.0L;
    SpectralReal target_scale = 0.0L;
    SpectralReal signed_bracket_ratio = 0.0L;
    SpectralReal first_square_ratio = 0.0L;
    SpectralReal absorption_ratio = 0.0L;
    bool finite = false;
};

class LocalSldRemainderAbsorptionObjective {
public:
    LocalSldRemainderAbsorptionObjective(
        const SpectralDynamics& dynamics,
        SpectralReal theta,
        TriadSelection selection =
            TriadSelection::local_without_equal_low_doubling());

    [[nodiscard]] LocalSldRemainderAbsorptionValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement absorption_ratio_gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    SpectralReal theta_ = 1.0L;
    TriadSelection selection_;
};

}  // namespace lemma
