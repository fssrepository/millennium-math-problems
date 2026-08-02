#pragma once

#include "projective_advection_decomposition.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveCoherenceObjectiveValue {
    std::size_t projective_shape_count = 0;
    SpectralReal coherent_norm2 = 0.0L;
    SpectralReal square_function_norm2 = 0.0L;
    SpectralReal synthesis_ratio = 0.0L;
    SpectralReal synthesis_amplification = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveCoherenceObjective {
public:
    LocalSldProjectiveCoherenceObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection);

    [[nodiscard]] LocalSldProjectiveCoherenceObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
};

}  // namespace lemma
