#pragma once

#include "projective_advection_decomposition.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveHeightStretchingObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    SpectralInteger shell_maximum_height = 0;
    std::size_t shell_shape_count = 0;
    SpectralReal signed_shell_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal aggregate_h1_norm2 = 0.0L;
    SpectralReal square_function_h1_norm2 = 0.0L;
    SpectralReal h1_synthesis_ratio = 0.0L;
    SpectralReal stretching_h1_alignment_squared = 0.0L;
    SpectralReal stretching_aware_h1_ratio = 0.0L;
    SpectralReal product_reconstruction_error = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveHeightStretchingObjective {
public:
    LocalSldProjectiveHeightStretchingObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveHeightStretchingObjectiveValue
    evaluate(const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    SpectralInteger core_maximum_height_ = 0;
    SpectralInteger shell_maximum_height_ = 0;
    int threads_ = 12;
};

}  // namespace lemma
