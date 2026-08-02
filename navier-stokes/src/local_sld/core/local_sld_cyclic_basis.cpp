#include "local_sld_cyclic_basis.hpp"

#include "spectral_galerkin.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

namespace lemma {

SpectralState LocalSldCyclicBasis::axis_state(int cutoff) {
    if (cutoff < 2) {
        throw std::invalid_argument(
            "the cyclic SLD basis requires cutoff at least two");
    }
    std::mt19937_64 generator(1);
    SpectralState state = SpectralStateFactory::random(cutoff, generator);
    for (ComplexVector& value : state.velocity) {
        value = {};
    }
    auto set_pair = [&](WaveVector wave, std::size_t component) {
        ComplexVector value{};
        value[component] = 1.0L;
        state.velocity[state.index.at(wave)] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    };
    set_pair(WaveVector{1, 0, 0}, 2);
    set_pair(WaveVector{0, 1, 0}, 0);
    set_pair(WaveVector{0, 0, 1}, 1);
    SpectralStateOps::normalize_energy(state);
    return state;
}

SpectralState LocalSldCyclicBasis::response_state(
    const SpectralDynamics& dynamics,
    const SpectralState& axis) {
    SpectralState response = axis;
    response.velocity = dynamics.advection_direct_partition(
        axis, TriadPartition::local);
    dynamics.enforce_constraints(response);
    SpectralStateOps::normalize_energy(response);
    return response;
}

SpectralState LocalSldCyclicBasis::cubic_response_state(
    const SpectralDynamics& dynamics,
    const SpectralState& axis,
    const SpectralState& response) {
    if (axis.waves != response.waves) {
        throw std::invalid_argument("cyclic basis layouts do not match");
    }
    SpectralState cubic = axis;
    cubic.velocity = dynamics.advection_jvp_direct(
        axis, response.velocity);
    dynamics.enforce_constraints(cubic);
    const SpectralReal axis_weight = pairing(
        axis.velocity, cubic.velocity);
    const SpectralReal response_weight = pairing(
        response.velocity, cubic.velocity);
    for (std::size_t mode = 0; mode < cubic.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            cubic.velocity[mode][component] -=
                axis_weight * axis.velocity[mode][component] +
                response_weight * response.velocity[mode][component];
        }
    }
    dynamics.enforce_constraints(cubic);
    SpectralStateOps::normalize_energy(cubic);
    return cubic;
}

SpectralState LocalSldCyclicBasis::mix(
    const SpectralState& axis,
    const SpectralState& response,
    SpectralReal angle) {
    if (axis.waves != response.waves) {
        throw std::invalid_argument("cyclic basis layouts do not match");
    }
    SpectralState state = axis;
    const SpectralReal left = std::cos(angle);
    const SpectralReal right = std::sin(angle);
    for (std::size_t mode = 0; mode < state.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            state.velocity[mode][component] =
                left * axis.velocity[mode][component] +
                right * response.velocity[mode][component];
        }
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

SpectralIncrement LocalSldCyclicBasis::angle_tangent(
    const SpectralState& axis,
    const SpectralState& response,
    SpectralReal angle) {
    if (axis.waves != response.waves) {
        throw std::invalid_argument("cyclic basis layouts do not match");
    }
    SpectralIncrement raw = axis.velocity;
    SpectralIncrement tangent = axis.velocity;
    const SpectralReal left = -std::sin(angle);
    const SpectralReal right = std::cos(angle);
    for (std::size_t mode = 0; mode < tangent.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            raw[mode][component] =
                std::cos(angle) * axis.velocity[mode][component] +
                std::sin(angle) * response.velocity[mode][component];
            tangent[mode][component] =
                left * axis.velocity[mode][component] +
                right * response.velocity[mode][component];
        }
    }
    const SpectralReal norm2 = pairing(raw, raw);
    if (!(norm2 > 0.0L)) {
        throw std::runtime_error("cyclic basis mixture has zero energy");
    }
    const SpectralReal inverse_norm = 1.0L / std::sqrt(norm2);
    const SpectralReal radial = pairing(raw, tangent) /
        (norm2 * std::sqrt(norm2));
    for (std::size_t mode = 0; mode < tangent.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            tangent[mode][component] =
                inverse_norm * tangent[mode][component] -
                radial * raw[mode][component];
        }
    }
    return tangent;
}

SpectralReal LocalSldCyclicBasis::pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("cyclic basis increments do not match");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

}  // namespace lemma
