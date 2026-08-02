#pragma once

#include "local_quartic_closure_objective.hpp"
#include "projective_advection_decomposition.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveCrossPowerObjectiveValue {
    std::size_t projective_shape_count = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal full_bracket = 0.0L;
    SpectralReal diagonal_bracket = 0.0L;
    SpectralReal cross_bracket = 0.0L;
    SpectralReal signed_cross_power_one = 0.0L;
    SpectralReal absolute_cross_power_one = 0.0L;
    SpectralReal squared_cross_power_one = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveCrossPowerObjective {
public:
    LocalSldProjectiveCrossPowerObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection);

    [[nodiscard]] LocalSldProjectiveCrossPowerObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

    [[nodiscard]] SpectralReal diagonal_bracket(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement diagonal_bracket_gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
};

}  // namespace lemma
