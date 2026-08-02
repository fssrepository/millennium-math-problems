#pragma once

#include "spectral_dynamics.hpp"
#include "triad_partition.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveNormalizationSchurObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    int row_shell = 0;
    std::size_t selected_shape_count = 0;
    std::size_t tail_shell_count = 0;
    SpectralReal full_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_aggregate_h1_norm2 = 0.0L;
    SpectralReal diagonal_tail_h2_norm2 = 0.0L;
    SpectralReal normalized_absolute_gram_row_sum = 0.0L;
    SpectralReal normalization_common_factor = 0.0L;
    SpectralReal schur_squared_majorant = 0.0L;
    SpectralReal height_half_compensated_schur_squared_majorant = 0.0L;
    bool finite = false;
};

// Smooth on each fixed Gram-sign stratum. It differentiates the exact
// fixed-row PNT-12 Schur quantity
// sqrt(H) * 9 S^2 B1 R_i D_H / (4 Z^3 P^5).
class LocalSldProjectiveNormalizationSchurObjective {
public:
    LocalSldProjectiveNormalizationSchurObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int row_shell,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveNormalizationSchurObjectiveValue
    evaluate(const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    SpectralInteger core_maximum_height_ = 0;
    int row_shell_ = 0;
    int threads_ = 12;
};

}  // namespace lemma
