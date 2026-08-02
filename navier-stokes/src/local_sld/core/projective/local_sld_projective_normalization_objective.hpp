#pragma once

#include "local_quartic_closure_objective.hpp"

namespace lemma {

struct LocalSldProjectiveNormalizationObjectiveValue {
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_stretching = 0.0L;
    SpectralReal selected_palinstrophy_cross = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal palinstrophy_normalization_power_one = 0.0L;
    SpectralReal squared_palinstrophy_normalization_power_one = 0.0L;
    bool finite = false;
};

// Targets the exact global term
// |S_full| |3 S_selected T_selected/(2P)| /(Z^2 P^2).
class LocalSldProjectiveNormalizationObjective {
public:
    LocalSldProjectiveNormalizationObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection);

    [[nodiscard]] LocalSldProjectiveNormalizationObjectiveValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
};

}  // namespace lemma
