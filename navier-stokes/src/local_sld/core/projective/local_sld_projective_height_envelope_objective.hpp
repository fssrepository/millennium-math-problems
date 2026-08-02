#pragma once

#include "local_quartic_closure_objective.hpp"

#include <array>
#include <cstddef>

namespace lemma {

struct LocalSldProjectiveHeightEnvelopeObjectiveValue {
    std::size_t height_shell_count = 0;
    std::size_t height_pair_count = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    std::array<SpectralReal, 5> absolute_component_brackets{};
    SpectralReal absolute_component_bracket_envelope = 0.0L;
    SpectralReal absolute_component_power_one_envelope = 0.0L;
    SpectralReal squared_component_power_one_envelope = 0.0L;
    bool finite = false;
    bool pairs_outer_and_advected_commutator = false;
};

class LocalSldProjectiveHeightEnvelopeObjective {
public:
    LocalSldProjectiveHeightEnvelopeObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        int threads = 12,
        bool pair_outer_and_advected_commutator = false);

    [[nodiscard]] LocalSldProjectiveHeightEnvelopeObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    int threads_ = 12;
    bool pair_outer_and_advected_commutator_ = false;
};

}  // namespace lemma
