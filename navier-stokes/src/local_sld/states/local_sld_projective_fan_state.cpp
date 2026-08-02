#include "local_sld_projective_fan_state.hpp"

#include "triad_partition.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>

namespace lemma {
namespace {

bool nonzero(const ComplexVector& value) {
    return std::norm(value[0]) + std::norm(value[1]) +
        std::norm(value[2]) > 0.0L;
}

void assign_reality_pair(
    SpectralState& state,
    WaveVector wave,
    ComplexVector value) {
    value = project_divergence_free(wave, value);
    state.velocity[state.index.at(wave)] = value;
    state.velocity[state.index.at(-wave)] = conjugate(value);
}

std::array<SpectralInteger, 3> primitive_signature(
    WaveVector p, WaveVector q, WaveVector target) {
    std::array<SpectralInteger, 3> signature{
        norm_squared(p), norm_squared(q), norm_squared(target)};
    std::sort(signature.begin(), signature.end());
    const SpectralInteger divisor = std::gcd(
        signature[0], std::gcd(signature[1], signature[2]));
    for (SpectralInteger& length : signature) {
        length /= divisor;
    }
    return signature;
}

}  // namespace

LocalSldProjectiveFanState LocalSldProjectiveFanStateFactory::make(
    int cutoff) {
    if (cutoff < 2 || cutoff > 48) {
        throw std::invalid_argument(
            "projective fan cutoff must be between 2 and 48");
    }
    LocalSldProjectiveFanState result;
    SpectralState& state = result.state;
    for (int z = -cutoff; z <= cutoff; ++z) {
        for (int y = -cutoff; y <= cutoff; ++y) {
            if (y != 0 || z != 0) {
                state.waves.push_back({0, y, z});
            }
        }
    }
    std::sort(state.waves.begin(), state.waves.end());
    state.velocity.resize(state.waves.size());
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        state.index.emplace(state.waves[index], index);
    }
    const WaveVector target{0, 0, cutoff};
    std::set<std::array<SpectralInteger, 3>> used_projective_shapes;
    for (int b = 1; 2 * b < cutoff; ++b) {
        for (int a = 1; a <= cutoff; ++a) {
            const WaveVector q{0, a, b};
            const WaveVector p{0, -a, cutoff - b};
            if (!TriadPartitioner::is_local(p, q, target)) {
                continue;
            }
            if (!used_projective_shapes.insert(
                     primitive_signature(p, q, target)).second) {
                continue;
            }
            if (nonzero(state.velocity[state.index.at(p)]) ||
                nonzero(state.velocity[state.index.at(q)])) {
                continue;
            }
            const SpectralReal p_norm = std::sqrt(
                static_cast<SpectralReal>(norm_squared(p)));
            assign_reality_pair(
                state, q,
                ComplexVector{
                    SpectralComplex{1.0L, 0.0L},
                    SpectralComplex{}, SpectralComplex{}});
            assign_reality_pair(
                state, p,
                ComplexVector{
                    SpectralComplex{},
                    SpectralComplex{
                        0.0L,
                        -static_cast<SpectralReal>(p.z) / p_norm},
                    SpectralComplex{
                        0.0L,
                        static_cast<SpectralReal>(p.y) / p_norm}});
            ++result.coherent_pairs;
        }
    }
    if (result.coherent_pairs == 0) {
        throw std::runtime_error(
            "projective fan cutoff produced no local coherent pairs");
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        if (is_positive_representative(state.waves[index]) &&
            nonzero(state.velocity[index])) {
            ++result.active_positive_modes;
        }
    }
    SpectralStateOps::normalize_energy(state);
    return result;
}

}  // namespace lemma
