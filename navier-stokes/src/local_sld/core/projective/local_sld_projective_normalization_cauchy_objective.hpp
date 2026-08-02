#pragma once

#include "projective_advection_decomposition.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveNormalizationCauchyObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    std::size_t selected_shape_count = 0;
    std::size_t tail_shape_count = 0;
    SpectralReal full_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_aggregate_h1_norm2 = 0.0L;
    SpectralReal tail_aggregate_h2_norm2 = 0.0L;
    SpectralReal cauchy_bound_power_one = 0.0L;
    SpectralReal squared_cauchy_bound_power_one = 0.0L;
    bool finite = false;
};

// Exact squared Cauchy majorant for the canonical selected-stretching x
// tail-palinstrophy-cross PNT channel:
// (9/4) S_full^2 ||A^(1/2)b_selected||_2^2
//       ||A b_tail||_2^2 / (Z^3 P^5).
class LocalSldProjectiveNormalizationCauchyObjective {
public:
    LocalSldProjectiveNormalizationCauchyObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveNormalizationCauchyObjectiveValue
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
