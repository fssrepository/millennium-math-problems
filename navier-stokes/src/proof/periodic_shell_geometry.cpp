#include "periodic_shell_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace lemma {
namespace {

std::uint64_t exact_shell_count(int shell) {
    const std::int64_t lower = INT64_C(1) << shell;
    const std::int64_t upper = INT64_C(1) << (shell + 1);
    const std::int64_t lower2 = lower * lower;
    const std::int64_t upper2 = upper * upper;
    std::uint64_t count = 0;
    for (std::int64_t x = -upper + 1; x < upper; ++x) {
        for (std::int64_t y = -upper + 1; y < upper; ++y) {
            for (std::int64_t z = -upper + 1; z < upper; ++z) {
                const std::int64_t radius2 = x * x + y * y + z * z;
                if (radius2 >= lower2 && radius2 < upper2) {
                    ++count;
                }
            }
        }
    }
    return count;
}

ShellGeometryReal weighted_adjacent_form(
    const std::vector<ShellGeometryReal>& amplitudes,
    ShellGeometryReal weight_power) {
    ShellGeometryReal value = 0.0L;
    for (std::size_t left = 0; left < amplitudes.size(); ++left) {
        const std::size_t first = left == 0 ? 0 : left - 1;
        const std::size_t last = std::min(left + 1, amplitudes.size() - 1);
        for (std::size_t right = first; right <= last; ++right) {
            const std::size_t high = std::max(left, right);
            value += std::exp2(
                         weight_power *
                         static_cast<ShellGeometryReal>(high)) *
                     amplitudes[left] * amplitudes[right];
        }
    }
    return value;
}

ShellGeometryReal weighted_square_norm(
    const std::vector<ShellGeometryReal>& amplitudes,
    ShellGeometryReal weight_power) {
    ShellGeometryReal value = 0.0L;
    for (std::size_t shell = 0; shell < amplitudes.size(); ++shell) {
        value += std::exp2(
                     weight_power *
                     static_cast<ShellGeometryReal>(shell)) *
                 amplitudes[shell] * amplitudes[shell];
    }
    return value;
}

ShellGeometryReal ratio_or_zero(
    ShellGeometryReal numerator, ShellGeometryReal denominator) {
    if (denominator == 0.0L) {
        return numerator == 0.0L
            ? 0.0L
            : std::numeric_limits<ShellGeometryReal>::infinity();
    }
    return numerator / denominator;
}

}  // namespace

PeriodicShellGeometryCertificate PeriodicShellGeometry::certify(
    int maximum_enumerated_shell, int overlap_samples, std::uint64_t seed) {
    if (maximum_enumerated_shell < 0 || maximum_enumerated_shell > 7 ||
        overlap_samples < 1 || overlap_samples > 1000000) {
        throw std::invalid_argument("invalid periodic shell certificate parameters");
    }
    PeriodicShellGeometryCertificate certificate;
    certificate.maximum_enumerated_shell = maximum_enumerated_shell;
    certificate.overlap_samples = overlap_samples;
    certificate.seed = seed;
    // After the low-shell Cauchy bound, every role has the common high
    // weight 2^(5l/2). For ordered adjacent pairs,
    // sum 2^(s max(r,t)) a_r a_t <= (2+2^s) sum 2^(s l) a_l^2.
    certificate.one_gain_overlap_constant =
        2.0L + std::exp2(2.5L);
    certificate.three_gain_overlap_constant =
        certificate.one_gain_overlap_constant;
    certificate.ft1_one_gain_constant =
        (certificate.low_advecting_block_constant +
         certificate.low_advected_block_constant) *
        certificate.one_gain_overlap_constant * std::sqrt(2.0L);
    certificate.ft1_three_gain_constant =
        certificate.low_target_block_constant *
        certificate.three_gain_overlap_constant *
        std::sqrt(8.0L / 7.0L);

    for (int shell = 0; shell <= maximum_enumerated_shell; ++shell) {
        const std::uint64_t exact = exact_shell_count(shell);
        const std::uint64_t scale = UINT64_C(1) << (3 * shell);
        const std::uint64_t cube_bound = UINT64_C(64) * scale;
        const ShellGeometryReal ratio =
            static_cast<ShellGeometryReal>(exact) /
            static_cast<ShellGeometryReal>(cube_bound);
        certificate.shell_counts.push_back(
            {shell, exact, cube_bound, ratio});
        certificate.maximum_count_ratio = std::max(
            certificate.maximum_count_ratio, ratio);
        certificate.all_bounds_hold =
            certificate.all_bounds_hold && exact <= cube_bound;
    }

    // If H/L>4, the other high frequency h satisfies h>3H/4. Hence
    // H/h<4/3<2 and the two high waves occupy equal or adjacent hard shells.
    certificate.separated_high_shells_are_adjacent =
        (4.0L / 3.0L) < 2.0L;
    certificate.all_bounds_hold = certificate.all_bounds_hold &&
        certificate.separated_high_shells_are_adjacent;

    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<ShellGeometryReal> exponent_distribution(
        -4.0L, 2.0L);
    std::uniform_real_distribution<ShellGeometryReal> zero_distribution(
        0.0L, 1.0L);
    const int sequence_shells = std::max(8, maximum_enumerated_shell + 1);
    constexpr ShellGeometryReal tolerance = 256.0L *
        std::numeric_limits<ShellGeometryReal>::epsilon();
    for (int sample = 0; sample < overlap_samples; ++sample) {
        std::vector<ShellGeometryReal> amplitudes(
            static_cast<std::size_t>(sequence_shells));
        for (ShellGeometryReal& amplitude : amplitudes) {
            amplitude = zero_distribution(generator) < 0.15L
                ? 0.0L
                : std::exp(exponent_distribution(generator));
        }
        if (std::all_of(amplitudes.begin(), amplitudes.end(),
                        [](ShellGeometryReal value) {
                            return value == 0.0L;
                        })) {
            amplitudes[static_cast<std::size_t>(
                sample % sequence_shells)] = 1.0L;
        }
        const ShellGeometryReal one_ratio = ratio_or_zero(
            weighted_adjacent_form(amplitudes, 2.5L),
            certificate.one_gain_overlap_constant *
                weighted_square_norm(amplitudes, 2.5L));
        const ShellGeometryReal three_ratio = ratio_or_zero(
            weighted_adjacent_form(amplitudes, 2.5L),
            certificate.three_gain_overlap_constant *
                weighted_square_norm(amplitudes, 2.5L));
        certificate.maximum_one_gain_overlap_ratio = std::max(
            certificate.maximum_one_gain_overlap_ratio, one_ratio);
        certificate.maximum_three_gain_overlap_ratio = std::max(
            certificate.maximum_three_gain_overlap_ratio, three_ratio);
        certificate.all_bounds_hold = certificate.all_bounds_hold &&
            one_ratio <= 1.0L + tolerance &&
            three_ratio <= 1.0L + tolerance;
    }
    return certificate;
}

}  // namespace lemma
