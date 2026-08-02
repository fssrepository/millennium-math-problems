#pragma once

#include "local_quartic_closure_objective.hpp"
#include "projective_core_family.hpp"

namespace lemma {

struct LocalSldProjectiveNormalizationObjectiveValue {
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_stretching = 0.0L;
    SpectralReal selected_palinstrophy_cross = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralInteger core_maximum_height = 0;
    std::size_t fixed_core_shape_count = 0;
    SpectralReal fixed_core_stretching = 0.0L;
    SpectralReal fixed_core_palinstrophy_cross = 0.0L;
    SpectralReal open_palinstrophy_normalization = 0.0L;
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
        TriadSelection selection,
        SpectralInteger core_maximum_height = 0,
        int threads = 1);

    [[nodiscard]] LocalSldProjectiveNormalizationObjectiveValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    SpectralInteger core_maximum_height_ = 0;
    int threads_ = 1;
    std::vector<ProjectivePrimitiveSignature> core_;
};

}  // namespace lemma
