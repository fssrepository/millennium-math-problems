#include "triad_partition.hpp"

#include <algorithm>
#include <limits>

namespace lemma {

bool TriadPartitioner::is_local(
    WaveVector first, WaveVector second, WaveVector third) {
    const SpectralInteger smallest = std::min({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    const SpectralInteger largest = std::max({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    return largest <= locality_ratio_squared * smallest;
}

int TriadPartitioner::dyadic_gap(
    WaveVector first, WaveVector second, WaveVector third) {
    const SpectralInteger smallest = std::min({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    const SpectralInteger largest = std::max({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    SpectralInteger upper = locality_ratio_squared * smallest;
    int gap = 0;
    while (largest > upper) {
        ++gap;
        if (upper > std::numeric_limits<SpectralInteger>::max() / 4) {
            break;
        }
        upper *= 4;
    }
    return gap;
}

bool TriadPartitioner::includes(
    WaveVector first, WaveVector second, WaveVector third,
    TriadPartition partition) {
    if (partition == TriadPartition::all) {
        return true;
    }
    const int gap = dyadic_gap(first, second, third);
    switch (partition) {
        case TriadPartition::all:
            return true;
        case TriadPartition::local:
            return gap == 0;
        case TriadPartition::nonlocal:
            return gap >= 1;
        case TriadPartition::near_nonlocal:
            return gap == 1;
        case TriadPartition::far_nonlocal:
            return gap >= 2;
    }
    return false;
}

bool TriadPartitioner::includes(
    const SpectralState& state, InteractionIndex interaction,
    TriadPartition partition) {
    const auto [first, second, third] = interaction;
    return includes(state.waves[first], state.waves[second],
                    state.waves[third], partition);
}

}  // namespace lemma
