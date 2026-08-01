#pragma once

#include "proof_scaling.hpp"

namespace lemma {

struct TransitionBlockScalingReport {
    Rational band_count_logarithm_power{1};
    Rational vortex_enstrophy_power{3, 4};
    Rational vortex_palinstrophy_power{3, 4};
    Rational young_remainder_conjugate{4};
    Rational post_young_logarithm_power{4};
    Rational post_young_enstrophy_power{3};
    Rational energy_time_integrable_enstrophy_power{1};
    Rational required_pointwise_depletion_power{1, 2};
    bool logarithmic_band_count_changes_polynomial_power = false;
    bool energy_identity_closes_transition_block = false;
};

class TransitionBlockScaling {
public:
    [[nodiscard]] static TransitionBlockScalingReport analyze();

    // Post-Young remainder divided by the energy-level Z term.
    [[nodiscard]] static long double normalized_remainder(long double z);
};

}  // namespace lemma
