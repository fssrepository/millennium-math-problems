#include "triad_partition.hpp"

#include <algorithm>

namespace lemma {

bool TriadPartitioner::is_local(
    WaveVector first, WaveVector second, WaveVector third) {
    const SpectralInteger smallest = std::min({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    const SpectralInteger largest = std::max({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    return largest <= locality_ratio_squared * smallest;
}

bool TriadPartitioner::includes(
    WaveVector first, WaveVector second, WaveVector third,
    TriadPartition partition) {
    if (partition == TriadPartition::all) {
        return true;
    }
    const bool local = is_local(first, second, third);
    return partition == TriadPartition::local ? local : !local;
}

bool TriadPartitioner::includes(
    const SpectralState& state, InteractionIndex interaction,
    TriadPartition partition) {
    const auto [first, second, third] = interaction;
    return includes(state.waves[first], state.waves[second],
                    state.waves[third], partition);
}

}  // namespace lemma
