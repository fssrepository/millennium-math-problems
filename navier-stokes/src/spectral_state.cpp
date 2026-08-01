#include "spectral_state.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace lemma {

WaveVector operator+(const WaveVector& a, const WaveVector& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

WaveVector operator-(const WaveVector& wave) {
    return {-wave.x, -wave.y, -wave.z};
}

SpectralInteger norm_squared(const WaveVector& wave) {
    return static_cast<SpectralInteger>(wave.x) * wave.x +
           static_cast<SpectralInteger>(wave.y) * wave.y +
           static_cast<SpectralInteger>(wave.z) * wave.z;
}

bool is_positive_representative(const WaveVector& wave) {
    if (wave.z != 0) {
        return wave.z > 0;
    }
    if (wave.y != 0) {
        return wave.y > 0;
    }
    return wave.x > 0;
}

SpectralComplex dot_hermitian(const ComplexVector& a, const ComplexVector& b) {
    SpectralComplex result{0.0L, 0.0L};
    for (std::size_t direction = 0; direction < 3; ++direction) {
        result += std::conj(a[direction]) * b[direction];
    }
    return result;
}

SpectralComplex wave_dot(const WaveVector& wave,
                         const ComplexVector& vector) {
    return static_cast<SpectralReal>(wave.x) * vector[0] +
           static_cast<SpectralReal>(wave.y) * vector[1] +
           static_cast<SpectralReal>(wave.z) * vector[2];
}

ComplexVector conjugate(const ComplexVector& vector) {
    return {std::conj(vector[0]), std::conj(vector[1]),
            std::conj(vector[2])};
}

ComplexVector project_divergence_free(const WaveVector& wave,
                                      ComplexVector vector) {
    const SpectralReal wave2 = static_cast<SpectralReal>(norm_squared(wave));
    const SpectralComplex coefficient = wave_dot(wave, vector) / wave2;
    vector[0] -= static_cast<SpectralReal>(wave.x) * coefficient;
    vector[1] -= static_cast<SpectralReal>(wave.y) * coefficient;
    vector[2] -= static_cast<SpectralReal>(wave.z) * coefficient;
    return vector;
}

int SpectralStateOps::cutoff(const SpectralState& state) {
    int result = 0;
    for (const WaveVector wave : state.waves) {
        result = std::max({result, std::abs(wave.x), std::abs(wave.y),
                           std::abs(wave.z)});
    }
    return result;
}

const std::vector<InteractionIndex>& SpectralStateOps::interactions(
    const SpectralState& state) {
    static std::map<int, std::vector<InteractionIndex>> cache;
    const int state_cutoff = cutoff(state);
    const auto existing = cache.find(state_cutoff);
    if (existing != cache.end()) {
        return existing->second;
    }
    std::vector<InteractionIndex> result;
    result.reserve(state.waves.size() * state.waves.size() / 2);
    for (std::size_t p = 0; p < state.waves.size(); ++p) {
        for (std::size_t q = 0; q < state.waves.size(); ++q) {
            const auto target = state.index.find(state.waves[p] + state.waves[q]);
            if (target != state.index.end()) {
                result.push_back({p, q, target->second});
            }
        }
    }
    return cache.emplace(state_cutoff, std::move(result)).first->second;
}

void SpectralStateOps::scale(SpectralState& state, SpectralReal factor) {
    for (auto& velocity : state.velocity) {
        for (SpectralComplex& component : velocity) {
            component *= factor;
        }
    }
}

SpectralReal SpectralStateOps::energy(const SpectralState& state) {
    SpectralReal result = 0.0L;
    for (const auto& velocity : state.velocity) {
        result += std::real(dot_hermitian(velocity, velocity));
    }
    return result;
}

void SpectralStateOps::normalize_energy(SpectralState& state,
                                        SpectralReal target_energy) {
    const SpectralReal current_energy = energy(state);
    if (!(current_energy > 0.0L)) {
        throw std::runtime_error("cannot normalize a zero spectral state");
    }
    scale(state, std::sqrt(target_energy / current_energy));
}

SpectralState SpectralStateFactory::random(int cutoff,
                                           std::mt19937_64& generator) {
    if (cutoff < 1 || cutoff > 12) {
        throw std::invalid_argument("--triad-cutoff must be between 1 and 12");
    }
    SpectralState state;
    for (int z = -cutoff; z <= cutoff; ++z) {
        for (int y = -cutoff; y <= cutoff; ++y) {
            for (int x = -cutoff; x <= cutoff; ++x) {
                if (x != 0 || y != 0 || z != 0) {
                    state.waves.push_back({x, y, z});
                }
            }
        }
    }
    std::sort(state.waves.begin(), state.waves.end());
    state.velocity.resize(state.waves.size());
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        state.index.emplace(state.waves[index], index);
    }
    std::normal_distribution<SpectralReal> normal(0.0L, 1.0L);
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        if (!is_positive_representative(wave)) {
            continue;
        }
        ComplexVector value{};
        for (SpectralComplex& component : value) {
            component = {normal(generator), normal(generator)};
        }
        value = project_divergence_free(wave, value);
        const SpectralReal decay =
            std::pow(1.0L + static_cast<SpectralReal>(norm_squared(wave)),
                     -1.25L);
        for (SpectralComplex& component : value) {
            component *= decay;
        }
        state.velocity[index] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    }
    return state;
}

namespace {

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

SpectralReal deterministic_signed_unit(std::uint64_t& state) {
    state = splitmix64(state);
    constexpr SpectralReal denominator =
        static_cast<SpectralReal>(1ULL << 53U);
    return 2.0L * static_cast<SpectralReal>(state >> 11U) / denominator - 1.0L;
}

std::uint64_t wave_seed(std::uint64_t seed, WaveVector wave) {
    std::uint64_t result = splitmix64(seed ^ 0xa0761d6478bd642fULL);
    result = splitmix64(result ^ static_cast<std::uint64_t>(
                                    static_cast<std::int64_t>(wave.x) + 0x10000));
    result = splitmix64(result ^ static_cast<std::uint64_t>(
                                    static_cast<std::int64_t>(wave.y) + 0x20000));
    return splitmix64(result ^ static_cast<std::uint64_t>(
                                  static_cast<std::int64_t>(wave.z) + 0x30000));
}

}  // namespace

SpectralState SpectralStateFactory::analytic(int cutoff, std::uint64_t seed,
                                             SpectralReal spectral_decay) {
    if (!(spectral_decay > 0.0L)) {
        throw std::invalid_argument("family spectral decay must be positive");
    }
    std::mt19937_64 layout_generator(0);
    SpectralState state = random(cutoff, layout_generator);
    for (auto& velocity : state.velocity) {
        velocity = {};
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        if (!is_positive_representative(wave)) {
            continue;
        }
        std::uint64_t mode_seed = wave_seed(seed, wave);
        ComplexVector value{};
        for (SpectralComplex& component : value) {
            component = {deterministic_signed_unit(mode_seed),
                         deterministic_signed_unit(mode_seed)};
        }
        value = project_divergence_free(wave, value);
        const SpectralReal amplitude = std::exp(
            -spectral_decay *
            std::sqrt(static_cast<SpectralReal>(norm_squared(wave))));
        for (SpectralComplex& component : value) {
            component *= amplitude;
        }
        state.velocity[index] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    }
    return state;
}

SpectralState SpectralStateFactory::mutate(const SpectralState& source,
                                           SpectralReal mutation,
                                           std::mt19937_64& generator,
                                           bool sparse) {
    SpectralState candidate = source;
    std::normal_distribution<SpectralReal> normal(0.0L, 1.0L);
    std::vector<std::size_t> positive_modes;
    for (std::size_t index = 0; index < candidate.waves.size(); ++index) {
        if (is_positive_representative(candidate.waves[index])) {
            positive_modes.push_back(index);
        }
    }
    std::vector<std::size_t> selected;
    if (sparse) {
        std::uniform_int_distribution<std::size_t> choose(0,
                                                          positive_modes.size() - 1);
        const std::size_t count = std::min<std::size_t>(
            positive_modes.size(), 1 + static_cast<std::size_t>(generator() % 3));
        while (selected.size() < count) {
            const std::size_t index = positive_modes[choose(generator)];
            if (std::find(selected.begin(), selected.end(), index) == selected.end()) {
                selected.push_back(index);
            }
        }
    } else {
        selected = positive_modes;
    }
    const SpectralReal per_mode =
        mutation / std::sqrt(static_cast<SpectralReal>(selected.size()));
    for (const std::size_t index : selected) {
        const WaveVector wave = candidate.waves[index];
        ComplexVector value = candidate.velocity[index];
        for (SpectralComplex& component : value) {
            component += per_mode *
                         SpectralComplex{normal(generator), normal(generator)};
        }
        value = project_divergence_free(wave, value);
        candidate.velocity[index] = value;
        candidate.velocity[candidate.index.at(-wave)] = conjugate(value);
    }
    SpectralStateOps::normalize_energy(candidate);
    return candidate;
}

SpectralState SpectralStateFactory::lift(const SpectralState& source,
                                         int target_cutoff,
                                         std::mt19937_64& generator) {
    SpectralState lifted = random(target_cutoff, generator);
    for (auto& velocity : lifted.velocity) {
        velocity = {};
    }
    for (std::size_t index = 0; index < source.waves.size(); ++index) {
        const auto target = lifted.index.find(source.waves[index]);
        if (target == lifted.index.end()) {
            throw std::runtime_error(
                "cannot lift a Fourier state to a smaller cutoff");
        }
        lifted.velocity[target->second] = source.velocity[index];
    }
    SpectralStateOps::normalize_energy(lifted);
    return lifted;
}

SpectralState SpectralStateFactory::project(const SpectralState& source,
                                            int target_cutoff) {
    if (target_cutoff < 1 ||
        target_cutoff >= SpectralStateOps::cutoff(source)) {
        throw std::invalid_argument(
            "spectral projection requires a smaller positive cutoff");
    }
    std::mt19937_64 layout_generator(0);
    SpectralState projected = random(target_cutoff, layout_generator);
    for (ComplexVector& value : projected.velocity) {
        value = {};
    }
    for (std::size_t index = 0; index < projected.waves.size(); ++index) {
        const auto source_index = source.index.find(projected.waves[index]);
        if (source_index != source.index.end()) {
            projected.velocity[index] =
                source.velocity[source_index->second];
        }
    }
    if (!(SpectralStateOps::energy(projected) > 0.0L)) {
        throw std::runtime_error(
            "spectral projection removed every nonzero mode");
    }
    SpectralStateOps::normalize_energy(projected);
    return projected;
}

}  // namespace lemma
