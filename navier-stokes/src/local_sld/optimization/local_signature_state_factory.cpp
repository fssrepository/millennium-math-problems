#include "local_signature_state_factory.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace lemma {

LocalSignatureStateProfile LocalSignatureStateFactory::parse(
    const std::string& name) {
    if (name == "decaying") {
        return LocalSignatureStateProfile::decaying;
    }
    if (name == "flat") {
        return LocalSignatureStateProfile::flat;
    }
    if (name == "outer-half-flat") {
        return LocalSignatureStateProfile::outer_half_flat;
    }
    throw std::invalid_argument(
        "signature profile must be decaying, flat, or outer-half-flat");
}

const char* LocalSignatureStateFactory::name(
    LocalSignatureStateProfile profile) {
    switch (profile) {
        case LocalSignatureStateProfile::decaying:
            return "decaying";
        case LocalSignatureStateProfile::flat:
            return "flat";
        case LocalSignatureStateProfile::outer_half_flat:
            return "outer-half-flat";
    }
    throw std::invalid_argument("invalid local signature state profile");
}

SpectralState LocalSignatureStateFactory::make(
    int cutoff, LocalSignatureStateProfile profile,
    std::uint64_t seed) {
    std::mt19937_64 generator(seed);
    SpectralState state = SpectralStateFactory::random(cutoff, generator);
    if (profile == LocalSignatureStateProfile::decaying) {
        SpectralStateOps::normalize_energy(state);
        return state;
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        if (profile == LocalSignatureStateProfile::outer_half_flat &&
            2 * std::max({std::abs(wave.x), std::abs(wave.y),
                          std::abs(wave.z)}) <= cutoff) {
            state.velocity[index] = {};
            continue;
        }
        const SpectralReal undo_decay = std::pow(
            1.0L + static_cast<SpectralReal>(norm_squared(wave)), 1.25L);
        for (SpectralComplex& component : state.velocity[index]) {
            component *= undo_decay;
        }
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

}  // namespace lemma
