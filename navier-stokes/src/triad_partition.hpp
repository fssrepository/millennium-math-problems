#pragma once

#include "spectral_state.hpp"

namespace lemma {

enum class TriadPartition { all, local, nonlocal };

class TriadPartitioner {
public:
    static constexpr SpectralInteger locality_ratio_squared = 4;

    [[nodiscard]] static bool is_local(
        WaveVector first, WaveVector second, WaveVector third);
    [[nodiscard]] static bool includes(
        WaveVector first, WaveVector second, WaveVector third,
        TriadPartition partition);
    [[nodiscard]] static bool includes(
        const SpectralState& state, InteractionIndex interaction,
        TriadPartition partition);
};

}  // namespace lemma
