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
    TriadSelection selection) {
    if (!selection.includes_gap(dyadic_gap(first, second, third))) {
        return false;
    }
    if (selection.signature_mode() ==
        TriadSelection::SignatureMode::none) {
        return true;
    }
    std::array<SpectralInteger, 3> signature{
        norm_squared(first), norm_squared(second), norm_squared(third)};
    std::sort(signature.begin(), signature.end());
    bool matches = false;
    if (selection.signature_mode() ==
            TriadSelection::SignatureMode::include_equal_low_doubling ||
        selection.signature_mode() ==
            TriadSelection::SignatureMode::exclude_equal_low_doubling) {
        matches = signature[0] == signature[1] &&
            signature[2] == 2 * signature[0];
        return selection.signature_mode() ==
                TriadSelection::SignatureMode::include_equal_low_doubling
            ? matches
            : !matches;
    }
    std::array<SpectralInteger, 3> requested =
        selection.squared_length_signature();
    std::sort(requested.begin(), requested.end());
    matches = signature == requested;
    return selection.signature_mode() ==
            TriadSelection::SignatureMode::include
        ? matches
        : !matches;
}

bool TriadPartitioner::includes(
    const SpectralState& state, InteractionIndex interaction,
    TriadSelection selection) {
    const auto [first, second, third] = interaction;
    return includes(state.waves[first], state.waves[second],
                    state.waves[third], selection);
}

}  // namespace lemma
