#include "proof_scaling.hpp"

#include <cstdlib>
#include <stdexcept>

namespace lemma {
namespace {

ScalingInteger gcd(ScalingInteger a, ScalingInteger b) {
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0) {
        const ScalingInteger remainder = a % b;
        a = b;
        b = remainder;
    }
    return a == 0 ? 1 : a;
}

}  // namespace

Rational::Rational(ScalingInteger numerator_in) : numerator(numerator_in) {}

Rational::Rational(ScalingInteger numerator_in, ScalingInteger denominator_in)
    : numerator(numerator_in), denominator(denominator_in) {
    if (denominator == 0) {
        throw std::invalid_argument("zero rational denominator");
    }
    if (denominator < 0) {
        numerator = -numerator;
        denominator = -denominator;
    }
    const ScalingInteger divisor = gcd(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;
}

long double Rational::value() const {
    return static_cast<long double>(numerator) /
           static_cast<long double>(denominator);
}

std::string Rational::str() const {
    if (denominator == 1) {
        return std::to_string(numerator);
    }
    return std::to_string(numerator) + "/" + std::to_string(denominator);
}

Rational operator+(const Rational& a, const Rational& b) {
    return {a.numerator * b.denominator + b.numerator * a.denominator,
            a.denominator * b.denominator};
}

Rational operator-(const Rational& a, const Rational& b) {
    return {a.numerator * b.denominator - b.numerator * a.denominator,
            a.denominator * b.denominator};
}

Rational operator*(const Rational& a, const Rational& b) {
    return {a.numerator * b.numerator, a.denominator * b.denominator};
}

Rational operator/(const Rational& a, const Rational& b) {
    if (b.numerator == 0) {
        throw std::invalid_argument("division by zero rational");
    }
    return {a.numerator * b.denominator, a.denominator * b.numerator};
}

bool operator==(const Rational& a, const Rational& b) {
    return a.numerator == b.numerator && a.denominator == b.denominator;
}

bool operator<(const Rational& a, const Rational& b) {
    return a.numerator * b.denominator < b.numerator * a.denominator;
}

bool operator<=(const Rational& a, const Rational& b) {
    return a < b || a == b;
}

ConcentrationScaling ScalingAnalyzer::analyze_concentration() {
    ConcentrationScaling result;
    const Rational amplitude_gain =
        Rational(2) * result.energy_normalizing_amplitude;
    result.fixed_energy_enstrophy = result.enstrophy + amplitude_gain;
    result.fixed_energy_pointwise_q =
        Rational(4) * result.depletion + result.fixed_energy_enstrophy;
    result.natural_critical_integrand =
        Rational(4) * result.depletion + Rational(2) * result.enstrophy;
    result.natural_integrated_l4 =
        result.natural_critical_integrand + result.time;
    result.pointwise_candidate_scale_compatible =
        result.fixed_energy_pointwise_q <= Rational(0);
    result.integrated_candidate_scale_critical =
        result.natural_integrated_l4 == Rational(0);
    return result;
}

StrongL4Reduction ScalingAnalyzer::analyze_strong_l4_reduction() {
    StrongL4Reduction result;
    result.exact_density_factorization =
        result.q_power_in_critical_density == Rational(1) &&
        result.enstrophy_power_in_critical_density == Rational(1) &&
        result.depletion_power_in_q == Rational(4) &&
        result.enstrophy_power_in_q == Rational(1);
    result.closes_integrated_l4_from_uniform_q =
        result.exact_density_factorization &&
        result.energy_bound_numerator == Rational(1) &&
        result.energy_bound_viscosity_denominator == Rational(2);
    return result;
}

ScalingCertificate ScalingAnalyzer::analyze_monomials(int denominator) {
    if (denominator < 4 || denominator > 100000) {
        throw std::invalid_argument(
            "--exponent-denominator must be between 4 and 100000");
    }
    ScalingCertificate result;
    const Rational zero(0);
    const Rational one(1);
    const Rational three_quarters(3, 4);
    for (int index = 0; index <= denominator; ++index) {
        const Rational a(index, denominator);
        const Rational b = three_quarters - Rational(2) * a;
        const Rational c = three_quarters + a;
        if (b < zero || c < zero || !(c < one)) {
            continue;
        }
        const Rational young_multiplier_power = one / (one - c);
        const Rational gamma = b * young_multiplier_power;
        const Rational pointwise_depletion =
            (gamma - one) / young_multiplier_power;
        const Rational energy_depletion =
            (gamma - Rational(2)) / young_multiplier_power;
        const ExponentCandidate candidate{a, b, c, young_multiplier_power, gamma,
                                          pointwise_depletion, energy_depletion};
        result.candidates.push_back(candidate);
        result.universal_quarter_depletion =
            result.universal_quarter_depletion &&
            energy_depletion == Rational(1, 4);
        if (!result.has_absorbable_candidate ||
            gamma < result.minimum_young_power) {
            result.minimum_young_power = gamma;
            result.minimizer = candidate;
        }
        result.has_absorbable_candidate = true;
        if (gamma <= one) {
            result.closing_candidate_exists = true;
        }
    }
    return result;
}

}  // namespace lemma
