#pragma once

#include "spectral_dynamics.hpp"

namespace lemma {

struct LocalQuarticClosureObjectiveValue {
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal signed_stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
    SpectralReal negative_commutator_pairing = 0.0L;
    SpectralReal advecting_slot = 0.0L;
    SpectralReal signed_two_entry_bracket = 0.0L;
    SpectralReal candidate_scale = 0.0L;
    SpectralReal constant_ratio = 0.0L;
    SpectralReal squared_constant_ratio = 0.0L;
    SpectralReal initial_frequency = 0.0L;
    SpectralReal initial_ep_shift = 0.0L;
    SpectralReal local_polynomial_numerator = 0.0L;
    SpectralReal local_polynomial_denominator = 0.0L;
    SpectralReal signed_local_sld_ratio = 0.0L;
    bool finite = false;
};

class LocalQuarticClosureObjective {
public:
    explicit LocalQuarticClosureObjective(
        const SpectralDynamics& dynamics);

    [[nodiscard]] LocalQuarticClosureObjectiveValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement two_entry_bracket_gradient(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement squared_constant_ratio_gradient(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement signed_local_sld_ratio_gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
};

}  // namespace lemma
