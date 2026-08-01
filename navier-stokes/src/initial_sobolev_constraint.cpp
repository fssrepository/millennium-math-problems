#include "initial_sobolev_constraint.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace lemma {

InitialSobolevConstraint::InitialSobolevConstraint(int order,
                                                   SpectralReal cap)
    : order_(order), cap_(cap) {
    if (order_ < 0 || order_ > 8 || !(cap_ >= 0.0L) ||
        !std::isfinite(cap_)) {
        throw std::invalid_argument("invalid initial Sobolev constraint");
    }
    if ((order_ == 0) != (cap_ == 0.0L)) {
        throw std::invalid_argument(
            "Sobolev order and cap must both be enabled or disabled");
    }
}

bool InitialSobolevConstraint::enabled() const {
    return order_ > 0;
}

int InitialSobolevConstraint::order() const {
    return order_;
}

SpectralReal InitialSobolevConstraint::cap() const {
    return cap_;
}

SpectralReal InitialSobolevConstraint::value(
    const SpectralState& state) const {
    if (!enabled()) {
        return 0.0L;
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        SpectralReal weight = 1.0L;
        for (int power = 0; power < order_; ++power) {
            weight *= wave2;
        }
        result += weight * std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
    }
    return result;
}

bool InitialSobolevConstraint::admissible(
    const SpectralState& state) const {
    return !enabled() || value(state) <= cap_ * (1.0L + 1e-12L);
}

SpectralIncrement InitialSobolevConstraint::energy_tangent_normal(
    const SpectralState& state) const {
    SpectralIncrement result(state.waves.size());
    if (!enabled()) {
        return result;
    }
    const SpectralReal energy = SpectralStateOps::energy(state);
    const SpectralReal moment = value(state);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        SpectralReal weight = 1.0L;
        for (int power = 0; power < order_; ++power) {
            weight *= wave2;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                (weight - moment / energy) *
                state.velocity[mode][component];
        }
    }
    return result;
}

void InitialSobolevConstraint::retract(
    SpectralState& state, SpectralReal target_energy) const {
    SpectralStateOps::normalize_energy(state, target_energy);
    if (!enabled() || admissible(state)) {
        return;
    }
    if (cap_ < target_energy * (1.0L - 1e-14L)) {
        throw std::invalid_argument(
            "Sobolev cap is below the minimum at the target energy");
    }
    const SpectralState unfiltered = state;
    auto filtered_value = [&](SpectralReal lambda, SpectralState* output) {
        SpectralState candidate = unfiltered;
        for (std::size_t mode = 0; mode < candidate.waves.size(); ++mode) {
            const SpectralReal wave2 = static_cast<SpectralReal>(
                norm_squared(candidate.waves[mode]));
            SpectralReal weight = 1.0L;
            for (int power = 0; power < order_; ++power) {
                weight *= wave2;
            }
            const SpectralReal factor =
                std::exp(-lambda * (weight - 1.0L));
            for (SpectralComplex& component : candidate.velocity[mode]) {
                component *= factor;
            }
        }
        SpectralStateOps::normalize_energy(candidate, target_energy);
        const SpectralReal result = value(candidate);
        if (output != nullptr) {
            *output = std::move(candidate);
        }
        return result;
    };

    SpectralReal lower = 0.0L;
    SpectralReal upper = 1e-12L;
    while (filtered_value(upper, nullptr) > cap_ && upper < 1e6L) {
        upper *= 2.0L;
    }
    if (filtered_value(upper, nullptr) > cap_) {
        throw std::runtime_error("cannot retract state to Sobolev cap");
    }
    for (int iteration = 0; iteration < 80; ++iteration) {
        const SpectralReal middle = 0.5L * (lower + upper);
        if (filtered_value(middle, nullptr) > cap_) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    static_cast<void>(filtered_value(upper, &state));
}

}  // namespace lemma
