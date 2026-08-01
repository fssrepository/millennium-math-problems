#pragma once

#include "spectral_state.hpp"

#include <limits>

namespace lemma {

enum class TriadPartition {
    all,
    local,
    nonlocal,
    near_nonlocal,
    far_nonlocal
};

class TriadSelection {
public:
    constexpr TriadSelection() = default;
    constexpr TriadSelection(TriadPartition partition) {
        switch (partition) {
            case TriadPartition::all:
                break;
            case TriadPartition::local:
                maximum_gap_ = 0;
                break;
            case TriadPartition::nonlocal:
                minimum_gap_ = 1;
                break;
            case TriadPartition::near_nonlocal:
                minimum_gap_ = 1;
                maximum_gap_ = 1;
                break;
            case TriadPartition::far_nonlocal:
                minimum_gap_ = 2;
                break;
        }
    }

    [[nodiscard]] static constexpr TriadSelection dyadic_tail(
        int minimum_gap) {
        return TriadSelection(minimum_gap, maximum_gap_value);
    }

    [[nodiscard]] constexpr int minimum_gap() const {
        return minimum_gap_;
    }
    [[nodiscard]] constexpr int maximum_gap() const {
        return maximum_gap_;
    }
    [[nodiscard]] constexpr bool is_all() const {
        return minimum_gap_ == 0 && maximum_gap_ == maximum_gap_value;
    }
    [[nodiscard]] constexpr bool includes_gap(int gap) const {
        return gap >= minimum_gap_ && gap <= maximum_gap_;
    }

private:
    static constexpr int maximum_gap_value =
        std::numeric_limits<int>::max();

    constexpr TriadSelection(int minimum_gap, int maximum_gap)
        : minimum_gap_(minimum_gap), maximum_gap_(maximum_gap) {}

    int minimum_gap_ = 0;
    int maximum_gap_ = maximum_gap_value;
};

class TriadPartitioner {
public:
    static constexpr SpectralInteger locality_ratio_squared = 4;

    [[nodiscard]] static bool is_local(
        WaveVector first, WaveVector second, WaveVector third);
    [[nodiscard]] static int dyadic_gap(
        WaveVector first, WaveVector second, WaveVector third);
    [[nodiscard]] static bool includes(
        WaveVector first, WaveVector second, WaveVector third,
        TriadSelection selection);
    [[nodiscard]] static bool includes(
        const SpectralState& state, InteractionIndex interaction,
        TriadSelection selection);
};

}  // namespace lemma
