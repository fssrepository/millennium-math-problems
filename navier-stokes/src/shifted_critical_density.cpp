#include "shifted_critical_density.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {

ShiftedCriticalDensityCertificate ShiftedCriticalDensityLemma::analyze() {
    const Rational energy_amplitude(2);
    const Rational enstrophy_amplitude(2);
    const Rational palinstrophy_amplitude(2);
    const Rational stretching_amplitude(3);

    const Rational energy_scaling(-1);
    const Rational enstrophy_scaling(1);
    const Rational palinstrophy_scaling(3);
    const Rational stretching_scaling(3);
    const Rational time_scaling(-2);

    ShiftedCriticalDensityCertificate result;
    result.density_amplitude_degree =
        Rational(4) * stretching_amplitude - enstrophy_amplitude -
        Rational(3) * palinstrophy_amplitude;
    result.initial_shift_amplitude_degree =
        energy_amplitude + palinstrophy_amplitude;
    result.density_scaling_exponent =
        Rational(4) * stretching_scaling - enstrophy_scaling -
        Rational(3) * palinstrophy_scaling;
    result.initial_shift_scaling_exponent =
        energy_scaling + palinstrophy_scaling;
    result.log_derivative_scaling_exponent =
        Rational(0) - time_scaling;
    const Rational initial_frequency_scaling =
        (enstrophy_scaling - energy_scaling) / Rational(2);
    result.gronwall_coefficient_scaling_exponent =
        initial_frequency_scaling + enstrophy_scaling;
    result.shift_matches_density =
        result.density_amplitude_degree ==
            result.initial_shift_amplitude_degree &&
        result.density_scaling_exponent ==
            result.initial_shift_scaling_exponent;
    result.gronwall_coefficient_is_critical =
        result.log_derivative_scaling_exponent ==
        result.gronwall_coefficient_scaling_exponent;
    result.energy_identity_closes_conditionally =
        result.shift_matches_density &&
        result.gronwall_coefficient_is_critical;
    return result;
}

ShiftedCriticalDensityDiagnostic ShiftedCriticalDensityAnalyzer::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralObjective& objective,
    const SpectralState& state,
    SpectralReal viscosity) {
    if (!(viscosity > 0.0L) || !std::isfinite(viscosity)) {
        throw std::invalid_argument(
            "shifted density diagnostic viscosity must be positive");
    }
    const StaticObjective local = objective.evaluate(
        state, TriadPartition::local);
    const SpectralIncrement gradient =
        objective.critical_integrand_gradient(
            state, TriadPartition::local);
    const SpectralIncrement velocity_derivative =
        dynamics.rhs(state, viscosity);
    ShiftedCriticalDensityDiagnostic result;
    result.energy = local.energy;
    result.enstrophy = local.enstrophy;
    result.palinstrophy = local.palinstrophy;
    result.local_critical_density = local.critical_integrand;
    result.initial_ep_shift = local.energy * local.palinstrophy;
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        result.density_derivative += std::real(
            dot_hermitian(gradient[mode], velocity_derivative[mode]));
    }
    const SpectralReal shifted_density =
        result.local_critical_density + result.initial_ep_shift;
    if (shifted_density > 0.0L && result.energy > 0.0L &&
        result.enstrophy > 0.0L) {
        result.shifted_log_derivative =
            result.density_derivative / shifted_density;
        result.initial_frequency = std::sqrt(
            result.enstrophy / result.energy);
        result.normalization =
            result.initial_frequency * result.enstrophy;
        if (result.normalization > 0.0L) {
            result.normalized_shifted_log_derivative =
                result.shifted_log_derivative / result.normalization;
        }
    }
    result.finite =
        std::isfinite(result.density_derivative) &&
        std::isfinite(result.shifted_log_derivative) &&
        std::isfinite(result.normalized_shifted_log_derivative);
    return result;
}

}  // namespace lemma
