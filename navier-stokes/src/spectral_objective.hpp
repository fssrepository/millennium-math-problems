#pragma once

#include "spectral_dynamics.hpp"

namespace lemma {

struct StaticObjective {
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal signed_vortex_stretching = 0.0L;
    SpectralReal vortex_stretching = 0.0L;
    SpectralReal depletion = 0.0L;
    SpectralReal energy_level_quantity = 0.0L;
    SpectralReal critical_integrand = 0.0L;
};

class SpectralObjective {
public:
    explicit SpectralObjective(const SpectralDynamics& dynamics);

    [[nodiscard]] StaticObjective evaluate(
        const SpectralState& state,
        TriadSelection selection = {}) const;
    [[nodiscard]] SpectralIncrement energy_level_gradient(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement critical_integrand_gradient(
        const SpectralState& state,
        TriadSelection selection = {}) const;

private:
    [[nodiscard]] SpectralIncrement signed_stretching_gradient(
        const SpectralState& state, TriadSelection selection) const;
    const SpectralDynamics& dynamics_;
};

}  // namespace lemma
