#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace lemma {

using ScalingInteger = std::int64_t;

class Rational {
public:
    ScalingInteger numerator = 0;
    ScalingInteger denominator = 1;

    Rational() = default;
    Rational(ScalingInteger numerator_in);
    Rational(ScalingInteger numerator_in, ScalingInteger denominator_in);

    [[nodiscard]] long double value() const;
    [[nodiscard]] std::string str() const;
};

Rational operator+(const Rational& a, const Rational& b);
Rational operator-(const Rational& a, const Rational& b);
Rational operator*(const Rational& a, const Rational& b);
Rational operator/(const Rational& a, const Rational& b);
bool operator==(const Rational& a, const Rational& b);
bool operator<(const Rational& a, const Rational& b);
bool operator<=(const Rational& a, const Rational& b);

struct ExponentCandidate {
    Rational energy;
    Rational enstrophy;
    Rational palinstrophy;
    Rational young_multiplier_power;
    Rational young_enstrophy_power;
    Rational pointwise_linear_depletion_power;
    Rational energy_integrable_depletion_power;
};

struct ScalingCertificate {
    std::vector<ExponentCandidate> candidates;
    Rational minimum_young_power{std::numeric_limits<ScalingInteger>::max(), 1};
    ExponentCandidate minimizer;
    bool has_absorbable_candidate = false;
    bool closing_candidate_exists = false;
    bool universal_quarter_depletion = true;
};

struct ConcentrationScaling {
    Rational energy{-1};
    Rational enstrophy{1};
    Rational palinstrophy{3};
    Rational vortex_stretching{3};
    Rational depletion{0};
    Rational time{-2};
    Rational energy_normalizing_amplitude{1, 2};
    Rational fixed_energy_enstrophy{2};
    Rational fixed_energy_pointwise_q{2};
    Rational natural_critical_integrand{2};
    Rational natural_integrated_l4{0};
    bool pointwise_candidate_scale_compatible = false;
    bool integrated_candidate_scale_critical = true;
};

struct StrongL4Reduction {
    Rational depletion_power_in_q{4};
    Rational enstrophy_power_in_q{1};
    Rational q_power_in_critical_density{1};
    Rational enstrophy_power_in_critical_density{1};
    Rational energy_bound_numerator{1};
    Rational energy_bound_viscosity_denominator{2};
    bool exact_density_factorization = true;
    bool closes_integrated_l4_from_uniform_q = true;
};

class ScalingAnalyzer {
public:
    static ScalingCertificate analyze_monomials(int denominator);
    static ConcentrationScaling analyze_concentration();
    static StrongL4Reduction analyze_strong_l4_reduction();
};

}  // namespace lemma
