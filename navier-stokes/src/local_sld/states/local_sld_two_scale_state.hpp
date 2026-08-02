#pragma once

#include "spectral_dynamics.hpp"

#include <vector>

namespace lemma {

class LocalSldTwoScaleState {
public:
    [[nodiscard]] static SpectralState cyclic_axis_mixture(
        int high_scale,
        SpectralReal high_to_low_energy_ratio);
    [[nodiscard]] static SpectralState cyclic_response_mixture(
        const SpectralDynamics& dynamics,
        int high_scale,
        SpectralReal high_to_low_energy_ratio,
        SpectralReal response_angle);
    [[nodiscard]] static std::vector<SpectralState>
    cyclic_response_mixtures(
        const SpectralDynamics& dynamics,
        int high_scale,
        SpectralReal high_to_low_energy_ratio,
        const std::vector<SpectralReal>& response_angles);
};

}  // namespace lemma
