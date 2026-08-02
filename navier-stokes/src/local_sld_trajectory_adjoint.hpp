#pragma once

#include "local_quartic_closure_objective.hpp"
#include "spectral_adjoint.hpp"

namespace lemma {

struct LocalSldTrajectoryValue {
    SpectralReal initial_frequency = 0.0L;
    SpectralReal initial_ep_shift = 0.0L;
    SpectralReal terminal_ratio = 0.0L;
    SpectralReal terminal_shift_fraction = 0.0L;
    int steps = 0;
    bool finite = false;
};

class LocalSldTrajectoryAdjoint {
public:
    explicit LocalSldTrajectoryAdjoint(
        const SpectralDynamics& dynamics,
        TriadSelection selection = TriadPartition::local);

    [[nodiscard]] LocalSldTrajectoryValue terminal_value(
        const SpectralState& initial,
        SpectralReal viscosity,
        SpectralReal time_step,
        int steps) const;
    [[nodiscard]] QTrajectoryGradient terminal_gradient(
        const SpectralState& initial,
        SpectralReal viscosity,
        SpectralReal time_step,
        int steps) const;
    [[nodiscard]] LocalSldTrajectoryValue maximum_value(
        const SpectralState& initial,
        SpectralReal viscosity,
        SpectralReal time_step,
        int steps) const;
    [[nodiscard]] QTrajectoryGradient maximum_gradient(
        const SpectralState& initial,
        SpectralReal viscosity,
        SpectralReal time_step,
        int steps) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
};

}  // namespace lemma
