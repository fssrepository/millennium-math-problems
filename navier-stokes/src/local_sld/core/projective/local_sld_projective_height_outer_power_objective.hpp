#pragma once

#include "local_quartic_closure_objective.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveHeightOuterPowerObjectiveValue {
    std::size_t height_shell_count = 0;
    std::size_t active_height_shell_count = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal diagonal_outer_h1_sum = 0.0L;
    SpectralReal signed_outer_power_one = 0.0L;
    SpectralReal absolute_outer_power_one = 0.0L;
    SpectralReal squared_outer_power_one = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveHeightOuterPowerObjective {
public:
    LocalSldProjectiveHeightOuterPowerObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveHeightOuterPowerObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    int threads_ = 12;
};

}  // namespace lemma
