#pragma once

#include "projective_advection_decomposition.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct ProjectiveQuarticDiagonalMoment {
    std::size_t projective_shape_count = 0;
    SpectralReal bracket = 0.0L;
    SpectralIncrement gradient;
};

class ProjectiveQuarticDiagonalKernel {
public:
    [[nodiscard]] static ProjectiveQuarticDiagonalMoment evaluate(
        const SpectralState& state,
        const std::vector<ProjectiveInteractionGroup>& groups,
        SpectralReal enstrophy,
        SpectralReal palinstrophy,
        bool compute_gradient);
};

}  // namespace lemma
