#include "far_tail_closure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace lemma {

FarTailClosureEvaluation FarTailClosure::evaluate(
    SpectralReal enstrophy, SpectralReal viscosity, int base_gap,
    const PeriodicShellGeometryCertificate& geometry) {
    if (!(viscosity > 0.0L) || !std::isfinite(viscosity) ||
        !geometry.all_bounds_hold) {
        throw std::invalid_argument("invalid far-tail closure parameters");
    }
    FarTailClosureEvaluation result;
    result.viscosity = viscosity;
    result.enstrophy = enstrophy;
    result.moving_gap = MovingGapController::decide(enstrophy, base_gap);
    const SpectralReal gap = static_cast<SpectralReal>(
        result.moving_gap.minimum_gap);
    result.tail_coefficient =
        geometry.ft1_one_gain_constant * std::exp2(-0.5L * gap) +
        geometry.ft1_three_gain_constant * std::exp2(-1.5L * gap);
    result.young_palinstrophy_coefficient = viscosity / 4.0L;
    result.young_remainder =
        27.0L / (4.0L * viscosity * viscosity * viscosity) *
        std::pow(result.tail_coefficient, 4.0L) *
        enstrophy * enstrophy * enstrophy;

    const SpectralReal base = static_cast<SpectralReal>(base_gap);
    result.linear_enstrophy_coefficient =
        54.0L / (viscosity * viscosity * viscosity) *
        (std::pow(geometry.ft1_one_gain_constant, 4.0L) *
             std::exp2(-2.0L * base) +
         std::pow(geometry.ft1_three_gain_constant, 4.0L) *
             std::exp2(-6.0L * base));
    result.linear_enstrophy_bound =
        result.linear_enstrophy_coefficient * enstrophy;
    if (result.linear_enstrophy_bound > 0.0L) {
        result.normalized_remainder_ratio =
            result.young_remainder / result.linear_enstrophy_bound;
    } else if (result.young_remainder > 0.0L) {
        result.normalized_remainder_ratio =
            std::numeric_limits<SpectralReal>::infinity();
    }
    return result;
}

FarTailClosureCertificate FarTailClosure::verify_random(
    SpectralReal viscosity, int base_gap, int samples, std::uint64_t seed,
    const PeriodicShellGeometryCertificate& geometry) {
    if (samples < 1 || samples > 1000000) {
        throw std::invalid_argument("invalid far-tail closure sample count");
    }
    FarTailClosureCertificate certificate;
    certificate.base_gap = base_gap;
    certificate.samples = samples;
    certificate.seed = seed;
    certificate.viscosity = viscosity;
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<SpectralReal> logarithm_distribution(
        -12.0L, 40.0L);
    constexpr SpectralReal tolerance = 512.0L *
        std::numeric_limits<SpectralReal>::epsilon();
    for (int sample = 0; sample < samples; ++sample) {
        const SpectralReal enstrophy = sample == 0
            ? 0.0L
            : std::exp2(logarithm_distribution(generator));
        const FarTailClosureEvaluation evaluation = evaluate(
            enstrophy, viscosity, base_gap, geometry);
        certificate.maximum_normalized_remainder_ratio = std::max(
            certificate.maximum_normalized_remainder_ratio,
            evaluation.normalized_remainder_ratio);
        certificate.all_bounds_hold = certificate.all_bounds_hold &&
            evaluation.normalized_remainder_ratio <= 1.0L + tolerance;
    }
    return certificate;
}

}  // namespace lemma
