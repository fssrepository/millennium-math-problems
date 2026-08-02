#pragma once

#include <array>
#include <complex>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <vector>

namespace lemma {

using SpectralInteger = std::int64_t;
using SpectralReal = long double;
using SpectralComplex = std::complex<SpectralReal>;

struct WaveVector {
    int x = 0;
    int y = 0;
    int z = 0;

    auto operator<=>(const WaveVector&) const = default;
};

WaveVector operator+(const WaveVector& a, const WaveVector& b);
WaveVector operator-(const WaveVector& wave);

using ComplexVector = std::array<SpectralComplex, 3>;
using SpectralIncrement = std::vector<ComplexVector>;

SpectralInteger norm_squared(const WaveVector& wave);
bool is_positive_representative(const WaveVector& wave);
SpectralComplex dot_hermitian(const ComplexVector& a, const ComplexVector& b);
SpectralComplex wave_dot(const WaveVector& wave, const ComplexVector& vector);
ComplexVector conjugate(const ComplexVector& vector);
ComplexVector project_divergence_free(const WaveVector& wave,
                                      ComplexVector vector);

struct SpectralState {
    std::vector<WaveVector> waves;
    std::vector<ComplexVector> velocity;
    std::map<WaveVector, std::size_t> index;
};

using InteractionIndex = std::array<std::size_t, 3>;

class SpectralStateOps {
public:
    static int cutoff(const SpectralState& state);
    static const std::vector<InteractionIndex>& interactions(
        const SpectralState& state);
    static void scale(SpectralState& state, SpectralReal factor);
    static SpectralReal energy(const SpectralState& state);
    static void normalize_energy(SpectralState& state,
                                 SpectralReal target_energy = 1.0L);
};

class SpectralStateFactory {
public:
    static SpectralState random(int cutoff, std::mt19937_64& generator);
    static SpectralState analytic(int cutoff, std::uint64_t seed,
                                  SpectralReal spectral_decay);
    static SpectralState mutate(const SpectralState& source,
                                SpectralReal mutation,
                                std::mt19937_64& generator, bool sparse);
    static SpectralState lift(const SpectralState& source, int target_cutoff,
                              std::mt19937_64& generator);
    static SpectralState project(const SpectralState& source,
                                 int target_cutoff);
};

}  // namespace lemma
