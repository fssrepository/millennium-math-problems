#pragma once

#include "spectral_dynamics.hpp"

namespace lemma {

struct LocalSldProjectedSquareReport {
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal hyperpalinstrophy = 0.0L;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
    SpectralReal expanded_negative_square = 0.0L;
    SpectralReal expanded_enstrophy_remainder = 0.0L;
    SpectralReal expanded_palinstrophy_remainder = 0.0L;
    SpectralReal expanded_total = 0.0L;
    SpectralReal completion_coefficient = 0.0L;
    SpectralReal completed_negative_square = 0.0L;
    SpectralReal completed_enstrophy_remainder = 0.0L;
    SpectralReal completed_hyperpalinstrophy_remainder = 0.0L;
    SpectralReal completed_total = 0.0L;
    SpectralReal completion_relative_error = 0.0L;
    SpectralReal absolute_target_scale_ratio = 0.0L;
    bool identity_verified = false;
};

// Exact finite-dimensional identity for any fixed triad selection:
// -<B,A B> + S^2/(2Z) + 3 S <A B,A u>/(2P)
// = -||A^(1/2)(B-c A u)||_2^2 + S^2/(2Z) + c^2 H3,
// where S=<A u,B>, c=3S/(4P), and H3=||A^(3/2)u||_2^2.
class LocalSldProjectedSquare {
public:
    [[nodiscard]] static LocalSldProjectedSquareReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        TriadSelection selection);
};

}  // namespace lemma
