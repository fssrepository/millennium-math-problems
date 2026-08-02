#pragma once

#include "local_quartic_closure_objective.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveHeightCommutatorRatioObjectiveValue {
    std::size_t height_shell_count = 0;
    SpectralReal commutator_paired_bracket_envelope = 0.0L;
    SpectralReal outer_h1_sum = 0.0L;
    SpectralReal coercivity_ratio = 0.0L;
    SpectralReal squared_coercivity_ratio = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveHeightCommutatorRatioObjective {
public:
    LocalSldProjectiveHeightCommutatorRatioObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        int threads = 12);

    [[nodiscard]]
    LocalSldProjectiveHeightCommutatorRatioObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    int threads_ = 12;
};

}  // namespace lemma
