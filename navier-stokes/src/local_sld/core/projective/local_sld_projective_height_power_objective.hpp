#pragma once

#include "local_quartic_closure_objective.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveHeightPowerObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    SpectralInteger shell_maximum_height = 0;
    std::size_t shell_shape_count = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal shell_internal_bracket = 0.0L;
    SpectralReal signed_shell_power_one = 0.0L;
    SpectralReal absolute_shell_power_one = 0.0L;
    SpectralReal squared_shell_power_one = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveHeightPowerObjective {
public:
    LocalSldProjectiveHeightPowerObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveHeightPowerObjectiveValue evaluate(
        const SpectralState& state) const;

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
