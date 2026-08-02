#include "projective_height_shell_partition.hpp"

#include <algorithm>
#include <stdexcept>

namespace lemma {

int ProjectiveHeightShellPartition::shell_index(
    SpectralInteger height) {
    if (height < 1) {
        throw std::invalid_argument(
            "projective primitive height must be positive");
    }
    int shell = 0;
    SpectralInteger upper = 1;
    while (upper < height) {
        upper *= 2;
        ++shell;
    }
    return shell;
}

SpectralInteger ProjectiveHeightShellPartition::minimum_height(
    int shell) {
    if (shell < 0) {
        throw std::invalid_argument(
            "projective height-shell index must be nonnegative");
    }
    if (shell == 0) {
        return 1;
    }
    SpectralInteger previous_upper = 1;
    for (int index = 1; index < shell; ++index) {
        previous_upper *= 2;
    }
    return previous_upper + 1;
}

SpectralInteger ProjectiveHeightShellPartition::maximum_height(
    int shell) {
    if (shell < 0) {
        throw std::invalid_argument(
            "projective height-shell index must be nonnegative");
    }
    SpectralInteger result = 1;
    for (int index = 0; index < shell; ++index) {
        result *= 2;
    }
    return result;
}

std::vector<ProjectiveHeightShellGroup>
ProjectiveHeightShellPartition::build(
    const std::vector<ProjectiveInteractionGroup>& groups) {
    int maximum_shell = 0;
    for (const auto& group : groups) {
        maximum_shell = std::max(
            maximum_shell,
            shell_index(group.primitive_squared_lengths[2]));
    }
    std::vector<ProjectiveHeightShellGroup> result(
        static_cast<std::size_t>(maximum_shell + 1));
    for (int shell = 0; shell <= maximum_shell; ++shell) {
        auto& target = result[static_cast<std::size_t>(shell)];
        target.shell = shell;
        target.minimum_height = minimum_height(shell);
        target.maximum_height = maximum_height(shell);
    }
    for (std::size_t index = 0; index < groups.size(); ++index) {
        auto& target = result[static_cast<std::size_t>(shell_index(
            groups[index].primitive_squared_lengths[2]))];
        target.group_indices.push_back(index);
        target.primitive_shapes.push_back(
            groups[index].primitive_squared_lengths);
        target.interaction_count += groups[index].interactions.size();
    }
    return result;
}

}  // namespace lemma
