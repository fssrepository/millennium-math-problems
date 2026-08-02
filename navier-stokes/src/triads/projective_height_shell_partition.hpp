#pragma once

#include "projective_advection_decomposition.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace lemma {

struct ProjectiveHeightShellGroup {
    int shell = 0;
    SpectralInteger minimum_height = 0;
    SpectralInteger maximum_height = 0;
    std::size_t interaction_count = 0;
    std::vector<std::size_t> group_indices;
    std::vector<std::array<SpectralInteger, 3>> primitive_shapes;
};

class ProjectiveHeightShellPartition {
public:
    [[nodiscard]] static int shell_index(SpectralInteger height);
    [[nodiscard]] static SpectralInteger minimum_height(int shell);
    [[nodiscard]] static SpectralInteger maximum_height(int shell);

    [[nodiscard]] static std::vector<ProjectiveHeightShellGroup> build(
        const std::vector<ProjectiveInteractionGroup>& groups);
};

}  // namespace lemma
