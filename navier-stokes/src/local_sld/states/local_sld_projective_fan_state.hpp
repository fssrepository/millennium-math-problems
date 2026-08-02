#pragma once

#include "spectral_state.hpp"

#include <cstddef>

namespace lemma {

struct LocalSldProjectiveFanState {
    SpectralState state;
    std::size_t coherent_pairs = 0;
    std::size_t active_positive_modes = 0;
};

class LocalSldProjectiveFanStateFactory {
public:
    [[nodiscard]] static LocalSldProjectiveFanState make(int cutoff);
};

}  // namespace lemma
