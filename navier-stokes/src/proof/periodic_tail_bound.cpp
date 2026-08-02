#include "periodic_tail_bound.hpp"

#include "triad_ledger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace lemma {

PeriodicTailBoundEvaluation PeriodicTailBound::evaluate(
    const SpectralState& state, int minimum_gap,
    const PeriodicShellGeometryCertificate& geometry) {
    if (minimum_gap < 2 || !geometry.all_bounds_hold) {
        throw std::invalid_argument("invalid periodic tail bound parameters");
    }
    PeriodicTailBoundEvaluation result;
    result.cutoff = SpectralStateOps::cutoff(state);
    result.minimum_gap = minimum_gap;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal amplitude2 = std::real(
            dot_hermitian(state.velocity[mode], state.velocity[mode]));
        result.enstrophy += wave2 * amplitude2;
        result.palinstrophy += wave2 * wave2 * amplitude2;
    }
    const TriadLedgerReport ledger = TriadLedger::analyze(state);
    for (const TriadGapLedgerRow& row : ledger.gaps) {
        if (row.dyadic_gap >= minimum_gap) {
            result.signed_tail_stretching += row.signed_stretching;
        }
    }
    result.absolute_signed_tail_stretching = std::abs(
        result.signed_tail_stretching);
    const SpectralReal common_norm =
        std::pow(result.enstrophy, 0.75L) *
        std::pow(result.palinstrophy, 0.75L);
    result.one_gain_bound =
        geometry.ft1_one_gain_constant *
        std::exp2(-0.5L * static_cast<SpectralReal>(minimum_gap)) *
        common_norm;
    result.three_gain_bound =
        geometry.ft1_three_gain_constant *
        std::exp2(-1.5L * static_cast<SpectralReal>(minimum_gap)) *
        common_norm;
    result.total_bound = result.one_gain_bound + result.three_gain_bound;
    if (result.total_bound > 0.0L) {
        result.bound_ratio =
            result.absolute_signed_tail_stretching / result.total_bound;
    } else if (result.absolute_signed_tail_stretching > 0.0L) {
        result.bound_ratio =
            std::numeric_limits<SpectralReal>::infinity();
    }
    return result;
}

PeriodicTailBoundCertificate PeriodicTailBound::verify_random(
    int cutoff, int minimum_gap, int samples, std::uint64_t seed,
    const PeriodicShellGeometryCertificate& geometry) {
    if (cutoff < 1 || cutoff > 16 || minimum_gap < 2 ||
        samples < 1 || samples > 100000) {
        throw std::invalid_argument(
            "invalid periodic tail random-certificate parameters");
    }
    PeriodicTailBoundCertificate certificate;
    certificate.cutoff = cutoff;
    certificate.minimum_gap = minimum_gap;
    certificate.samples = samples;
    certificate.seed = seed;
    std::mt19937_64 generator(seed);
    constexpr SpectralReal tolerance = 256.0L *
        std::numeric_limits<SpectralReal>::epsilon();
    for (int sample = 0; sample < samples; ++sample) {
        SpectralState state = SpectralStateFactory::random(cutoff, generator);
        SpectralStateOps::normalize_energy(state);
        const PeriodicTailBoundEvaluation evaluation = evaluate(
            state, minimum_gap, geometry);
        certificate.maximum_bound_ratio = std::max(
            certificate.maximum_bound_ratio, evaluation.bound_ratio);
        certificate.nonzero_tail_seen = certificate.nonzero_tail_seen ||
            evaluation.absolute_signed_tail_stretching > 0.0L;
        certificate.all_bounds_hold = certificate.all_bounds_hold &&
            evaluation.bound_ratio <= 1.0L + tolerance;
    }
    certificate.all_bounds_hold = certificate.all_bounds_hold &&
        certificate.nonzero_tail_seen;
    return certificate;
}

}  // namespace lemma
