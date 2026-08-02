#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <string>

namespace lemma {

enum class LocalSignatureStateProfile {
    decaying,
    flat,
    outer_half_flat
};

class LocalSignatureStateFactory {
public:
    [[nodiscard]] static LocalSignatureStateProfile parse(
        const std::string& name);
    [[nodiscard]] static const char* name(
        LocalSignatureStateProfile profile);
    [[nodiscard]] static SpectralState make(
        int cutoff, LocalSignatureStateProfile profile,
        std::uint64_t seed);
};

}  // namespace lemma
