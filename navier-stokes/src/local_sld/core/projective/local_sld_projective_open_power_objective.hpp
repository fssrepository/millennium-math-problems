#pragma once

#include "local_quartic_closure_objective.hpp"
#include "projective_core_family.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct LocalSldProjectiveOpenPowerObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    std::size_t fixed_core_shape_count = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal selected_bracket = 0.0L;
    SpectralReal fixed_core_bracket = 0.0L;
    SpectralReal open_bracket = 0.0L;
    SpectralReal signed_open_power_one = 0.0L;
    SpectralReal absolute_open_power_one = 0.0L;
    SpectralReal squared_open_power_one = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveOpenPowerObjective {
public:
    LocalSldProjectiveOpenPowerObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int threads = 1);

    [[nodiscard]] LocalSldProjectiveOpenPowerObjectiveValue evaluate(
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
