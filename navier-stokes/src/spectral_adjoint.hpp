#pragma once

#include "spectral_objective.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct QTrajectoryGradient {
    SpectralIncrement initial_gradient;
    SpectralReal objective_value = 0.0L;
    int objective_step = 0;
    int total_steps = 0;
    std::size_t checkpoint_count = 0;
};

class SpectralAdjoint {
public:
    SpectralAdjoint(const SpectralDynamics& dynamics,
                    const SpectralObjective& objective);

    [[nodiscard]] QTrajectoryGradient terminal_q_gradient(
        const SpectralState& initial, SpectralReal viscosity,
        SpectralReal time_step, int steps) const;
    [[nodiscard]] QTrajectoryGradient maximum_q_gradient(
        const SpectralState& initial, SpectralReal viscosity,
        SpectralReal time_step, int steps) const;
    [[nodiscard]] QTrajectoryGradient q_gain_gradient(
        const SpectralState& initial, SpectralReal viscosity,
        SpectralReal time_step, int steps) const;
    [[nodiscard]] QTrajectoryGradient q_increase_gradient(
        const SpectralState& initial, SpectralReal viscosity,
        SpectralReal time_step, int steps) const;

private:
    [[nodiscard]] QTrajectoryGradient reverse_from_step(
        const std::vector<SpectralState>& checkpoints,
        SpectralReal viscosity, SpectralReal time_step,
        int objective_step) const;

    const SpectralDynamics& dynamics_;
    const SpectralObjective& objective_;
};

}  // namespace lemma
