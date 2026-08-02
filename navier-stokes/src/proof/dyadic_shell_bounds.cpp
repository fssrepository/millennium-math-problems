#include "dyadic_shell_bounds.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace lemma {
namespace {

ShellReal safe_ratio(ShellReal numerator, ShellReal denominator) {
    if (denominator == 0.0L) {
        return numerator == 0.0L
            ? 0.0L
            : std::numeric_limits<ShellReal>::infinity();
    }
    return numerator / denominator;
}

ShellReal dyadic_power(ShellReal exponent) {
    return std::exp2(exponent);
}

}  // namespace

DyadicShellBoundReport DyadicShellBounds::analyze(
    const std::vector<ShellReal>& amplitudes, int minimum_gap) {
    if (amplitudes.empty() || minimum_gap < 1) {
        throw std::invalid_argument(
            "dyadic shell amplitudes must be nonempty and gap positive");
    }
    DyadicShellBoundReport report;
    report.shells = static_cast<int>(amplitudes.size());
    report.minimum_gap = minimum_gap;
    report.low_one_derivative_constant = std::sqrt(2.0L);
    report.low_three_derivative_constant = std::sqrt(8.0L / 7.0L);

    for (std::size_t shell = 0; shell < amplitudes.size(); ++shell) {
        const ShellReal amplitude = amplitudes[shell];
        if (!(amplitude >= 0.0L) || !std::isfinite(amplitude)) {
            throw std::invalid_argument(
                "dyadic shell amplitudes must be finite and nonnegative");
        }
        const ShellReal j = static_cast<ShellReal>(shell);
        const ShellReal amplitude2 = amplitude * amplitude;
        report.enstrophy += dyadic_power(2.0L * j) * amplitude2;
        report.palinstrophy += dyadic_power(4.0L * j) * amplitude2;
        report.high_moment += dyadic_power(2.5L * j) * amplitude2;
    }
    report.high_moment_bound =
        std::pow(report.enstrophy, 0.75L) *
        std::pow(report.palinstrophy, 0.25L);
    report.high_moment_ratio = safe_ratio(
        report.high_moment, report.high_moment_bound);

    ShellReal low_one_sum = 0.0L;
    ShellReal low_three_sum = 0.0L;
    for (std::size_t high_shell = 0;
         high_shell < amplitudes.size(); ++high_shell) {
        if (high_shell >= static_cast<std::size_t>(minimum_gap)) {
            const std::size_t added_shell =
                high_shell - static_cast<std::size_t>(minimum_gap);
            const ShellReal low = static_cast<ShellReal>(added_shell);
            low_one_sum += dyadic_power(2.5L * low) *
                           amplitudes[added_shell];
            low_three_sum += dyadic_power(3.5L * low) *
                             amplitudes[added_shell];
        }

        const ShellReal high = static_cast<ShellReal>(high_shell);
        const ShellReal sqrt_palinstrophy =
            std::sqrt(report.palinstrophy);
        const ShellReal low_one_bound =
            report.low_one_derivative_constant *
            dyadic_power(0.5L *
                         (high - static_cast<ShellReal>(minimum_gap))) *
            sqrt_palinstrophy;
        const ShellReal low_three_bound =
            report.low_three_derivative_constant *
            dyadic_power(1.5L *
                         (high - static_cast<ShellReal>(minimum_gap))) *
            sqrt_palinstrophy;
        report.maximum_low_one_derivative_ratio = std::max(
            report.maximum_low_one_derivative_ratio,
            safe_ratio(low_one_sum, low_one_bound));
        report.maximum_low_three_derivative_ratio = std::max(
            report.maximum_low_three_derivative_ratio,
            safe_ratio(low_three_sum, low_three_bound));

        const ShellReal high_amplitude2 =
            amplitudes[high_shell] * amplitudes[high_shell];
        report.one_gain_tail +=
            dyadic_power(2.0L * high) * high_amplitude2 * low_one_sum;
        report.three_gain_tail +=
            dyadic_power(high) * high_amplitude2 * low_three_sum;
    }

    const ShellReal common_norm =
        std::pow(report.enstrophy, 0.75L) *
        std::pow(report.palinstrophy, 0.75L);
    report.one_gain_tail_bound =
        report.low_one_derivative_constant *
        dyadic_power(-0.5L * static_cast<ShellReal>(minimum_gap)) *
        common_norm;
    report.three_gain_tail_bound =
        report.low_three_derivative_constant *
        dyadic_power(-1.5L * static_cast<ShellReal>(minimum_gap)) *
        common_norm;
    report.one_gain_tail_ratio = safe_ratio(
        report.one_gain_tail, report.one_gain_tail_bound);
    report.three_gain_tail_ratio = safe_ratio(
        report.three_gain_tail, report.three_gain_tail_bound);
    return report;
}

DyadicShellRandomCertificate DyadicShellBounds::verify_random(
    int shells, int minimum_gap, int samples, std::uint64_t seed) {
    if (shells < 1 || shells > 1024 || minimum_gap < 1 ||
        minimum_gap > shells || samples < 1 || samples > 1000000) {
        throw std::invalid_argument("invalid dyadic shell verification parameters");
    }
    DyadicShellRandomCertificate certificate;
    certificate.shells = shells;
    certificate.minimum_gap = minimum_gap;
    certificate.samples = samples;
    certificate.seed = seed;
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<ShellReal> exponent_distribution(
        -3.0L, 1.0L);
    std::uniform_real_distribution<ShellReal> zero_distribution(0.0L, 1.0L);
    constexpr ShellReal tolerance = 256.0L *
        std::numeric_limits<ShellReal>::epsilon();

    for (int sample = 0; sample < samples; ++sample) {
        std::vector<ShellReal> amplitudes(static_cast<std::size_t>(shells));
        for (ShellReal& amplitude : amplitudes) {
            amplitude = zero_distribution(generator) < 0.15L
                ? 0.0L
                : std::exp(exponent_distribution(generator));
        }
        if (std::all_of(amplitudes.begin(), amplitudes.end(),
                        [](ShellReal value) { return value == 0.0L; })) {
            amplitudes[static_cast<std::size_t>(sample % shells)] = 1.0L;
        }
        const DyadicShellBoundReport report = analyze(
            amplitudes, minimum_gap);
        certificate.maximum_high_moment_ratio = std::max(
            certificate.maximum_high_moment_ratio,
            report.high_moment_ratio);
        certificate.maximum_low_one_derivative_ratio = std::max(
            certificate.maximum_low_one_derivative_ratio,
            report.maximum_low_one_derivative_ratio);
        certificate.maximum_low_three_derivative_ratio = std::max(
            certificate.maximum_low_three_derivative_ratio,
            report.maximum_low_three_derivative_ratio);
        certificate.maximum_one_gain_tail_ratio = std::max(
            certificate.maximum_one_gain_tail_ratio,
            report.one_gain_tail_ratio);
        certificate.maximum_three_gain_tail_ratio = std::max(
            certificate.maximum_three_gain_tail_ratio,
            report.three_gain_tail_ratio);
        certificate.all_bounds_hold = certificate.all_bounds_hold &&
            report.high_moment_ratio <= 1.0L + tolerance &&
            report.maximum_low_one_derivative_ratio <= 1.0L + tolerance &&
            report.maximum_low_three_derivative_ratio <= 1.0L + tolerance &&
            report.one_gain_tail_ratio <= 1.0L + tolerance &&
            report.three_gain_tail_ratio <= 1.0L + tolerance;
    }
    return certificate;
}

}  // namespace lemma
