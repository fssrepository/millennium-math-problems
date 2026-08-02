#pragma once

#include "projective_advection_decomposition.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveNormalizationAlignmentObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    std::size_t selected_shape_count = 0;
    std::size_t tail_shape_count = 0;
    SpectralReal selected_stretching = 0.0L;
    SpectralReal tail_palinstrophy_cross = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_aggregate_h1_norm2 = 0.0L;
    SpectralReal tail_aggregate_h2_norm2 = 0.0L;
    SpectralReal selected_stretching_h1_alignment_squared = 0.0L;
    SpectralReal tail_palinstrophy_cross_h2_alignment_squared = 0.0L;
    SpectralReal normalization_alignment_product_squared = 0.0L;
    bool finite = false;
};

// Exact squared actual/Cauchy ratio for the canonical
// selected-stretching x tail-palinstrophy-cross channel:
// s^2 t_tail^2 / (Z P ||A^(1/2)b||_2^2 ||A b_tail||_2^2).
class LocalSldProjectiveNormalizationAlignmentObjective {
public:
    LocalSldProjectiveNormalizationAlignmentObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveNormalizationAlignmentObjectiveValue
    evaluate(const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    SpectralInteger core_maximum_height_ = 0;
    int threads_ = 12;
};

}  // namespace lemma
