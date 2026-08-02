#include "local_sld_two_scale_state.hpp"

#include "local_sld_cyclic_basis.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralState dilated_mixture(
    const SpectralState& low,
    int high_scale,
    SpectralReal high_to_low_energy_ratio) {
    SpectralState high = low;
    for (ComplexVector& value : high.velocity) {
        value = {};
    }
    for (std::size_t mode = 0; mode < low.waves.size(); ++mode) {
        const ComplexVector value = low.velocity[mode];
        if (!(std::real(dot_hermitian(value, value)) > 0.0L)) {
            continue;
        }
        const WaveVector source = low.waves[mode];
        const WaveVector target{
            high_scale * source.x,
            high_scale * source.y,
            high_scale * source.z};
        const auto found = high.index.find(target);
        if (found == high.index.end()) {
            throw std::invalid_argument(
                "dilated cyclic response exceeds the requested cutoff");
        }
        high.velocity[found->second] = value;
    }
    SpectralStateOps::normalize_energy(high);
    SpectralState result = low;
    const SpectralReal low_weight = 1.0L / std::sqrt(
        1.0L + high_to_low_energy_ratio);
    const SpectralReal high_weight = std::sqrt(
        high_to_low_energy_ratio /
        (1.0L + high_to_low_energy_ratio));
    for (std::size_t mode = 0; mode < result.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            result.velocity[mode][component] =
                low_weight * low.velocity[mode][component] +
                high_weight * high.velocity[mode][component];
        }
    }
    SpectralStateOps::normalize_energy(result);
    return result;
}

}  // namespace

SpectralState LocalSldTwoScaleState::cyclic_axis_mixture(
    int high_scale,
    SpectralReal high_to_low_energy_ratio) {
    if (high_scale < 2 || high_scale > 12 ||
        !(high_to_low_energy_ratio > 0.0L) ||
        !std::isfinite(high_to_low_energy_ratio)) {
        throw std::invalid_argument(
            "two-scale cyclic state requires scale 2..12 and positive finite energy ratio");
    }
    const SpectralState low =
        LocalSldCyclicBasis::axis_state(high_scale);
    return dilated_mixture(
        low, high_scale, high_to_low_energy_ratio);
}

SpectralState LocalSldTwoScaleState::cyclic_response_mixture(
    const SpectralDynamics& dynamics,
    int high_scale,
    SpectralReal high_to_low_energy_ratio,
    SpectralReal response_angle) {
    if (high_scale < 2 || high_scale > 12 ||
        !(high_to_low_energy_ratio > 0.0L) ||
        !std::isfinite(high_to_low_energy_ratio) ||
        !std::isfinite(response_angle)) {
        throw std::invalid_argument(
            "two-scale cyclic response requires valid scale, energy ratio, and angle");
    }
    const SpectralState axis =
        LocalSldCyclicBasis::axis_state(high_scale);
    const SpectralState response =
        LocalSldCyclicBasis::response_state(dynamics, axis);
    const SpectralState low = LocalSldCyclicBasis::mix(
        axis, response, response_angle);
    return dilated_mixture(
        low, high_scale, high_to_low_energy_ratio);
}

std::vector<SpectralState>
LocalSldTwoScaleState::cyclic_response_mixtures(
    const SpectralDynamics& dynamics,
    int high_scale,
    SpectralReal high_to_low_energy_ratio,
    const std::vector<SpectralReal>& response_angles) {
    if (high_scale < 2 || high_scale > 12 ||
        !(high_to_low_energy_ratio > 0.0L) ||
        !std::isfinite(high_to_low_energy_ratio)) {
        throw std::invalid_argument(
            "two-scale cyclic response family requires valid scale and energy ratio");
    }
    const SpectralState axis =
        LocalSldCyclicBasis::axis_state(high_scale);
    const SpectralState response =
        LocalSldCyclicBasis::response_state(dynamics, axis);
    std::vector<SpectralState> result;
    result.reserve(response_angles.size());
    for (const SpectralReal angle : response_angles) {
        if (!std::isfinite(angle)) {
            throw std::invalid_argument(
                "two-scale response family angle must be finite");
        }
        const SpectralState low = LocalSldCyclicBasis::mix(
            axis, response, angle);
        result.push_back(dilated_mixture(
            low, high_scale, high_to_low_energy_ratio));
    }
    return result;
}

}  // namespace lemma
