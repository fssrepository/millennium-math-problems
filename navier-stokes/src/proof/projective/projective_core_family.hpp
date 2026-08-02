#pragma once

#include "spectral_state.hpp"

#include <array>
#include <vector>

namespace lemma {

using ProjectivePrimitiveSignature =
    std::array<SpectralInteger, 3>;

class ProjectiveCoreFamily {
public:
    [[nodiscard]] static std::vector<ProjectivePrimitiveSignature>
    through_maximum_height(SpectralInteger maximum_height);
};

}  // namespace lemma
