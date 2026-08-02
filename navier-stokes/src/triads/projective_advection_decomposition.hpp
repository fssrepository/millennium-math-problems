#pragma once

#include "spectral_state.hpp"
#include "triad_partition.hpp"

#include <array>
#include <vector>

namespace lemma {

struct ProjectiveInteractionGroup {
    std::array<SpectralInteger, 3> primitive_squared_lengths{};
    std::vector<InteractionIndex> interactions;
};

class ProjectiveAdvectionDecomposition {
public:
    [[nodiscard]] static std::vector<ProjectiveInteractionGroup> group(
        const SpectralState& state,
        TriadSelection selection);

    [[nodiscard]] static SpectralIncrement evaluate(
        const SpectralState& state,
        const ProjectiveInteractionGroup& group);

    [[nodiscard]] static SpectralIncrement vjp(
        const SpectralState& state,
        const ProjectiveInteractionGroup& group,
        const SpectralIncrement& output_cotangent);
};

}  // namespace lemma
