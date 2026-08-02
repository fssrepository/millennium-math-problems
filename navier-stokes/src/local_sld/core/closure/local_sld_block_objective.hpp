#pragma once

#include "local_quartic_closure_objective.hpp"

namespace lemma {

enum class LocalSldBlock {
    full,
    selected_closed,
    complement_closed,
    mixed
};

struct LocalSldBlockObjectiveValue {
    SpectralReal full_constant_ratio = 0.0L;
    SpectralReal selected_constant_ratio = 0.0L;
    SpectralReal complement_constant_ratio = 0.0L;
    SpectralReal mixed_constant_ratio = 0.0L;
    SpectralReal block_constant_ratio = 0.0L;
    SpectralReal normalized_stretching = 0.0L;
    SpectralReal common_shape_factor = 0.0L;
    SpectralReal block_sld_ratio = 0.0L;
    SpectralReal full_sld_ratio = 0.0L;
    SpectralReal ratio_reconstruction_error = 0.0L;
    bool finite = false;
};

class LocalSldBlockObjective {
public:
    LocalSldBlockObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selected,
        LocalSldBlock block);

    [[nodiscard]] LocalSldBlockObjectiveValue evaluate(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selected_;
    TriadSelection complement_;
    LocalSldBlock block_;
};

}  // namespace lemma
