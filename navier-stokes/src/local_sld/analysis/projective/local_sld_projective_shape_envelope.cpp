#include "local_sld_projective_shape_envelope.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace lemma {
namespace {

int dyadic_level(SpectralInteger value) {
    if (value <= 0) {
        throw std::invalid_argument(
            "projective primitive height must be positive");
    }
    int level = 0;
    while (value >= 2) {
        value /= 2;
        ++level;
    }
    return level;
}

SpectralReal relative_error(SpectralReal left, SpectralReal right) {
    return std::abs(left - right) /
        std::max({std::abs(left), std::abs(right), 1e-30L});
}

SpectralReal fitted_shell_slope(
    const std::vector<LocalSldProjectiveShapeShellRow>& rows) {
    SpectralReal count = 0.0L;
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    for (const auto& row : rows) {
        if (!(row.absolute_power_one_sum > 0.0L)) {
            continue;
        }
        const SpectralReal x = static_cast<SpectralReal>(
            row.dyadic_height_level);
        const SpectralReal y = std::log2(row.absolute_power_one_sum);
        count += 1.0L;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }
    const SpectralReal denominator = count * sum_xx - sum_x * sum_x;
    if (count < 2.0L || std::abs(denominator) < 1e-30L) {
        return 0.0L;
    }
    return (count * sum_xy - sum_x * sum_y) / denominator;
}

struct ShellAccumulator {
    LocalSldProjectiveShapeShellRow row;
    SpectralReal weighted_sine_sum = 0.0L;
    SpectralReal weighted_aspect_sum = 0.0L;
};

}  // namespace

LocalSldProjectiveShapeEnvelopeReport
LocalSldProjectiveShapeEnvelope::analyze(
    const LocalSldRemainderProjectiveReport& projective) {
    LocalSldProjectiveShapeEnvelopeReport report;
    report.projective_shape_count = projective.projective_shape_count;
    report.expected_signed_total =
        projective.reconstructed_power_one_total;
    std::map<int, ShellAccumulator> shell_map;
    SpectralReal weighted_sine_sum = 0.0L;
    SpectralReal weighted_aspect_sum = 0.0L;
    SpectralReal weighted_span_sum = 0.0L;
    for (const LocalSldRemainderProjectiveEntry& shape :
         projective.shapes) {
        const SpectralInteger a = shape.primitive_squared_lengths[0];
        const SpectralInteger b = shape.primitive_squared_lengths[1];
        const SpectralInteger c = shape.primitive_squared_lengths[2];
        if (a <= 0 || b <= 0 || c <= 0) {
            throw std::invalid_argument(
                "projective shape envelope requires positive lengths");
        }
        const SpectralReal ar = static_cast<SpectralReal>(a);
        const SpectralReal br = static_cast<SpectralReal>(b);
        const SpectralReal cr = static_cast<SpectralReal>(c);
        const SpectralReal cosine = (cr - ar - br) /
            (2.0L * std::sqrt(ar * br));
        const SpectralReal sine_squared = std::max(
            0.0L, 1.0L - cosine * cosine);
        const SpectralReal aspect = std::sqrt(cr / ar);
        const SpectralReal dyadic_span = std::log2(aspect);
        const SpectralReal weight = std::abs(shape.power_one_total);
        const int level = dyadic_level(c);
        ShellAccumulator& shell = shell_map[level];
        shell.row.dyadic_height_level = level;
        shell.row.minimum_primitive_height =
            shell.row.minimum_primitive_height == 0
                ? c
                : std::min(shell.row.minimum_primitive_height, c);
        shell.row.maximum_primitive_height = std::max(
            shell.row.maximum_primitive_height, c);
        ++shell.row.shape_count;
        shell.row.signed_power_one_total += shape.power_one_total;
        shell.row.absolute_power_one_sum += weight;
        shell.row.squared_power_one_sum +=
            shape.power_one_total * shape.power_one_total;
        shell.weighted_sine_sum += weight * sine_squared;
        shell.weighted_aspect_sum += weight * aspect;
        report.reconstructed_signed_total += shape.power_one_total;
        report.absolute_total += weight;
        weighted_sine_sum += weight * sine_squared;
        weighted_aspect_sum += weight * aspect;
        weighted_span_sum += weight * dyadic_span;
        report.primitive_height_half_moment +=
            weight * std::sqrt(cr);
        report.primitive_height_first_moment += weight * cr;
    }
    report.height_shells.reserve(shell_map.size());
    for (auto& [level, shell] : shell_map) {
        static_cast<void>(level);
        if (shell.row.squared_power_one_sum > 0.0L) {
            shell.row.effective_shapes =
                shell.row.absolute_power_one_sum *
                shell.row.absolute_power_one_sum /
                shell.row.squared_power_one_sum;
        }
        if (shell.row.absolute_power_one_sum > 0.0L) {
            shell.row.signed_alignment =
                std::abs(shell.row.signed_power_one_total) /
                shell.row.absolute_power_one_sum;
            shell.row.weighted_angle_sine_squared =
                shell.weighted_sine_sum /
                shell.row.absolute_power_one_sum;
            shell.row.weighted_length_aspect_ratio =
                shell.weighted_aspect_sum /
                shell.row.absolute_power_one_sum;
        }
        if (report.absolute_total > 0.0L) {
            shell.row.absolute_fraction =
                shell.row.absolute_power_one_sum / report.absolute_total;
        }
        report.height_shells.push_back(shell.row);
    }
    if (report.absolute_total > 0.0L) {
        report.contribution_weighted_angle_sine_squared =
            weighted_sine_sum / report.absolute_total;
        report.contribution_weighted_length_aspect_ratio =
            weighted_aspect_sum / report.absolute_total;
        report.contribution_weighted_dyadic_span =
            weighted_span_sum / report.absolute_total;
    }
    report.fitted_absolute_height_shell_slope = fitted_shell_slope(
        report.height_shells);
    report.reconstruction_error = relative_error(
        report.reconstructed_signed_total,
        report.expected_signed_total);
    report.exact_reconstruction = projective.exact_reconstruction &&
        report.reconstruction_error < 1e-13L;
    return report;
}

}  // namespace lemma
