#include "local_sld_cyclic_orbit_basis.hpp"

#include <array>
#include <cmath>
#include <random>
#include <stdexcept>

namespace lemma {
namespace {

int sign(int value) {
    return value < 0 ? -1 : 1;
}

ComplexVector transverse_value(WaveVector wave) {
    const int absolute_x = std::abs(wave.x);
    const int absolute_y = std::abs(wave.y);
    const int absolute_z = std::abs(wave.z);
    std::array<int, 3> rotated_wave{};
    int major = 0;
    if (absolute_x == 2 && absolute_y == 1 && absolute_z == 1) {
        rotated_wave = {wave.x, wave.y, wave.z};
        major = 0;
    } else if (absolute_y == 2 && absolute_z == 1 && absolute_x == 1) {
        rotated_wave = {wave.y, wave.z, wave.x};
        major = 1;
    } else if (absolute_z == 2 && absolute_x == 1 && absolute_y == 1) {
        rotated_wave = {wave.z, wave.x, wave.y};
        major = 2;
    } else {
        return {};
    }
    constexpr SpectralReal a = 0.0L;
    constexpr SpectralReal b = 1.0L;
    constexpr SpectralReal c = -1.0L;
    const std::array<SpectralReal, 3> rotated_value{
        -a * static_cast<SpectralReal>(sign(rotated_wave[1])),
        b * static_cast<SpectralReal>(sign(rotated_wave[0])),
        c * static_cast<SpectralReal>(
            sign(rotated_wave[0]) * sign(rotated_wave[1]) *
            sign(rotated_wave[2]))};
    std::array<SpectralReal, 3> value{};
    if (major == 0) {
        value = rotated_value;
    } else if (major == 1) {
        value = {
            rotated_value[2], rotated_value[0], rotated_value[1]};
    } else {
        value = {
            rotated_value[1], rotated_value[2], rotated_value[0]};
    }
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    return {
        imaginary_unit * value[0],
        imaginary_unit * value[1],
        imaginary_unit * value[2]};
}

}  // namespace

SpectralState LocalSldCyclicOrbitBasis::transverse_two_one_one(
    int cutoff) {
    if (cutoff < 2 || cutoff > 12) {
        throw std::invalid_argument(
            "the transverse (2,1,1) orbit requires cutoff 2..12");
    }
    std::mt19937_64 generator(211);
    SpectralState state = SpectralStateFactory::random(cutoff, generator);
    for (ComplexVector& value : state.velocity) {
        value = {};
    }
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const WaveVector wave = state.waves[mode];
        if (!is_positive_representative(wave)) {
            continue;
        }
        ComplexVector value = transverse_value(wave);
        value = project_divergence_free(wave, value);
        state.velocity[mode] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

namespace {

SpectralState three_one_zero(int cutoff, bool forward) {
    if (cutoff < 3 || cutoff > 12) {
        throw std::invalid_argument(
            "the (3,1,0) cyclic orbits require cutoff 3..12");
    }
    std::mt19937_64 generator(forward ? 310 : 301);
    SpectralState state = SpectralStateFactory::random(cutoff, generator);
    for (ComplexVector& value : state.velocity) {
        value = {};
    }
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const WaveVector wave = state.waves[mode];
        if (!is_positive_representative(wave)) {
            continue;
        }
        const std::array<int, 3> components{wave.x, wave.y, wave.z};
        int major = -1;
        for (int component = 0; component < 3; ++component) {
            if (std::abs(components[static_cast<std::size_t>(component)]) ==
                3) {
                major = component;
            }
        }
        if (major < 0) {
            continue;
        }
        const int minor = forward ? (major + 1) % 3 : (major + 2) % 3;
        const int zero = forward ? (major + 2) % 3 : (major + 1) % 3;
        if (std::abs(components[static_cast<std::size_t>(minor)]) != 1 ||
            components[static_cast<std::size_t>(zero)] != 0) {
            continue;
        }
        ComplexVector value{};
        const SpectralReal amplitude = forward
            ? static_cast<SpectralReal>(
                  sign(components[static_cast<std::size_t>(major)]))
            : -static_cast<SpectralReal>(
                  sign(components[static_cast<std::size_t>(minor)]));
        value[static_cast<std::size_t>(zero)] =
            imaginary_unit * amplitude;
        value = project_divergence_free(wave, value);
        state.velocity[mode] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

}  // namespace

SpectralState LocalSldCyclicOrbitBasis::forward_three_one_zero(
    int cutoff) {
    return three_one_zero(cutoff, true);
}

SpectralState LocalSldCyclicOrbitBasis::backward_three_one_zero(
    int cutoff) {
    return three_one_zero(cutoff, false);
}

}  // namespace lemma
