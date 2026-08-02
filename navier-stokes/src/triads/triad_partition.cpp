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
            signature[2] == 2 * signature[0] &&
            signature[0] >= selection.minimum_low_squared() &&
            signature[0] <
                selection.maximum_low_squared_exclusive();
        return selection.signature_mode() ==
                TriadSelection::SignatureMode::include_equal_low_doubling
            ? matches
            : !matches;
    }
    std::array<SpectralInteger, 3> requested =
        selection.squared_length_signature();
    std::sort(requested.begin(), requested.end());
    matches = signature == requested;
    if (selection.signature_mode() ==
        TriadSelection::SignatureMode::
            exclude_equal_low_doubling_and_signature) {
        const bool equal_low_doubling =
            signature[0] == signature[1] &&
            signature[2] == 2 * signature[0];
        return !equal_low_doubling && !matches;
    }
    if (selection.signature_mode() ==
            TriadSelection::SignatureMode::
                include_equal_low_double_triple ||
        selection.signature_mode() ==
            TriadSelection::SignatureMode::
                exclude_equal_low_double_triple ||
        selection.signature_mode() ==
            TriadSelection::SignatureMode::
                exclude_equal_low_double_triple_and_signature) {
        const bool double_or_triple =
            signature[0] == signature[1] &&
            (signature[2] == 2 * signature[0] ||
             signature[2] == 3 * signature[0]);
        if (selection.signature_mode() ==
            TriadSelection::SignatureMode::
                include_equal_low_double_triple) {
            return double_or_triple;
        }
        return !double_or_triple &&
            (selection.signature_mode() !=
                 TriadSelection::SignatureMode::
                     exclude_equal_low_double_triple_and_signature ||
             !matches);
    }
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
