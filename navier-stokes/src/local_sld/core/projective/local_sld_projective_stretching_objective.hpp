#pragma once

#include "projective_advection_decomposition.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveStretchingObjectiveValue {
    std::size_t projective_shape_count = 0;
    SpectralReal signed_selected_stretching = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal coherent_norm2 = 0.0L;
    SpectralReal square_function_norm2 = 0.0L;
    SpectralReal coherent_synthesis_ratio = 0.0L;
    SpectralReal stretching_alignment_squared = 0.0L;
    SpectralReal stretching_aware_synthesis_ratio = 0.0L;
    SpectralReal product_reconstruction_error = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveStretchingObjective {
public:
    LocalSldProjectiveStretchingObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection);

    [[nodiscard]] LocalSldProjectiveStretchingObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
};

}  // namespace lemma
