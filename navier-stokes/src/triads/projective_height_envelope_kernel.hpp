#pragma once

#include "projective_advection_decomposition.hpp"
#include "triad_partition.hpp"

#include <array>
#include <cstddef>

namespace lemma {

struct ProjectiveHeightEnvelopeMoment {
    std::size_t height_shell_count = 0;
    std::size_t height_pair_count = 0;
    std::array<SpectralReal, 5> absolute_component_sums{};
    SpectralReal absolute_component_envelope = 0.0L;
    SpectralIncrement gradient;
};

// Exact value and a sign-chamber gradient of the five-component absolute
// envelope of the complete dyadic primitive-height quartet matrix. At a zero
// component the selected subgradient is zero.
class ProjectiveHeightEnvelopeKernel {
public:
    [[nodiscard]] static ProjectiveHeightEnvelopeMoment evaluate(
        const SpectralState& state,
        TriadSelection selection,
        SpectralReal enstrophy,
        SpectralReal palinstrophy,
        bool compute_gradient,
        int threads = 1,
        bool pair_outer_and_advected_commutator = false);
};

}  // namespace lemma
