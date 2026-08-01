#pragma once

#include "helical_sector_objective.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

struct HelicalSectorTrajectoryGradient {
    SpectralIncrement initial_gradient;
    SpectralReal objective_value = 0.0L;
    int total_steps = 0;
    std::size_t checkpoint_count = 0;
};

class HelicalSectorAdjoint {
public:
    explicit HelicalSectorAdjoint(const SpectralDynamics& dynamics)
        : dynamics_(dynamics) {}

    [[nodiscard]] HelicalSectorTrajectoryGradient critical_integral_gradient(
        const SpectralState& initial, SpectralReal viscosity,
        SpectralReal time_step, int steps,
        HelicalSectorSelection selection) const;
    [[nodiscard]] SpectralReal critical_integral(
        const SpectralState& initial, SpectralReal viscosity,
        SpectralReal time_step, int steps,
        HelicalSectorSelection selection) const;

private:
    const SpectralDynamics& dynamics_;
};

}  // namespace lemma
