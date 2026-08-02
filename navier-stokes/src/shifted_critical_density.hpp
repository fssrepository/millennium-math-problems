#pragma once

#include "proof_scaling.hpp"
#include "spectral_objective.hpp"

namespace lemma {

struct ShiftedCriticalDensityCertificate {
    Rational density_amplitude_degree;
    Rational initial_shift_amplitude_degree;
    Rational density_scaling_exponent;
    Rational initial_shift_scaling_exponent;
    Rational log_derivative_scaling_exponent;
    Rational gronwall_coefficient_scaling_exponent;
    bool shift_matches_density = false;
    bool gronwall_coefficient_is_critical = false;
    bool energy_identity_closes_conditionally = false;
};

class ShiftedCriticalDensityLemma {
public:
    // B0=E(0)P(0). If the still-open differential inequality
    // d/dt log(C_local+B0) <= A k0 Z holds uniformly in the cutoff, then
    // the energy identity closes the local critical integral.
    [[nodiscard]] static ShiftedCriticalDensityCertificate analyze();
};

struct ShiftedCriticalDensityDiagnostic {
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal local_critical_density = 0.0L;
    SpectralReal initial_ep_shift = 0.0L;
    SpectralReal density_derivative = 0.0L;
    SpectralReal shifted_log_derivative = 0.0L;
    SpectralReal initial_frequency = 0.0L;
    SpectralReal normalization = 0.0L;
    SpectralReal normalized_shifted_log_derivative = 0.0L;
    bool finite = false;
};

class ShiftedCriticalDensityAnalyzer {
public:
    [[nodiscard]] static ShiftedCriticalDensityDiagnostic evaluate(
        const SpectralDynamics& dynamics,
        const SpectralObjective& objective,
        const SpectralState& state,
        SpectralReal viscosity);
};

}  // namespace lemma
