#pragma once

#include "local_quartic_closure_objective.hpp"

namespace lemma {

struct LocalSldShapePowerObjectiveValue {
    int power = 0;
    SpectralReal bracket_constant_ratio = 0.0L;
    SpectralReal normalized_stretching = 0.0L;
    SpectralReal absolute_power_product = 0.0L;
    SpectralReal squared_power_product = 0.0L;
    bool finite = false;
};

class LocalSldShapePowerObjective {
public:
    LocalSldShapePowerObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        int power);

    [[nodiscard]] LocalSldShapePowerObjectiveValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    int power_ = 0;
};

}  // namespace lemma
