#include "local_sld_projective_height_tail_summary.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal relative_error(
    SpectralReal computed,
    SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

}  // namespace

LocalSldProjectiveHeightTailReport
LocalSldProjectiveHeightTailSummary::summarize(
    const LocalSldProjectiveHeightMatrixReport& matrix) {
    if (matrix.shells.empty()) {
        throw std::invalid_argument(
            "projective height-tail summary requires height shells");
    }
    LocalSldProjectiveHeightTailReport report;
    for (std::size_t core = 0; core < matrix.shells.size(); ++core) {
        LocalSldProjectiveHeightTailRow row;
        row.last_core_shell = static_cast<int>(core);
        row.core_maximum_height = matrix.shells[core].maximum_height;
        for (std::size_t shell = 0; shell < matrix.shells.size(); ++shell) {
            if (shell <= core) {
                row.core_stretching += matrix.shells[shell].stretching;
                row.core_palinstrophy_cross +=
                    matrix.shells[shell].palinstrophy_cross;
            } else {
                row.tail_stretching += matrix.shells[shell].stretching;
                row.tail_palinstrophy_cross +=
                    matrix.shells[shell].palinstrophy_cross;
            }
        }
        if (matrix.selected_palinstrophy > 0.0L) {
            const SpectralReal factor =
                1.5L * matrix.power_one_scale /
                matrix.selected_palinstrophy;
            row.core_stretching_tail_cross_power_one = factor *
                row.core_stretching * row.tail_palinstrophy_cross;
            row.tail_stretching_core_cross_power_one = factor *
                row.tail_stretching * row.core_palinstrophy_cross;
            row.tail_stretching_tail_cross_power_one = factor *
                row.tail_stretching * row.tail_palinstrophy_cross;
        }
        SpectralReal open_square_sum = 0.0L;
        for (const auto& entry : matrix.entries) {
            const bool first_core =
                entry.first_shell <= row.last_core_shell;
            const bool second_core =
                entry.second_shell <= row.last_core_shell;
            if (first_core && second_core) {
                row.core_aggregate_h1_norm2 +=
                    entry.aggregate_h1_pairing;
                row.core_aggregate_h2_norm2 +=
                    entry.aggregate_h2_pairing;
                row.core_internal_power_one += entry.power_one;
            } else if (first_core) {
                row.core_tail_power_one += entry.power_one;
                row.open_absolute_power_one_sum +=
                    std::abs(entry.power_one);
                open_square_sum += entry.power_one * entry.power_one;
            } else {
                row.tail_aggregate_h1_norm2 +=
                    entry.aggregate_h1_pairing;
                row.tail_aggregate_h2_norm2 +=
                    entry.aggregate_h2_pairing;
                row.tail_internal_power_one += entry.power_one;
                row.open_absolute_power_one_sum +=
                    std::abs(entry.power_one);
                open_square_sum += entry.power_one * entry.power_one;
            }
            if (!second_core) {
                row.open_outer_square_power_one +=
                    entry.outer_square * matrix.power_one_scale;
                row.open_advected_commutator_power_one +=
                    entry.advected_commutator * matrix.power_one_scale;
                row.open_advecting_nested_power_one +=
                    entry.advecting_nested * matrix.power_one_scale;
                row.open_enstrophy_normalization_power_one +=
                    entry.enstrophy_normalization *
                    matrix.power_one_scale;
                row.open_palinstrophy_normalization_power_one +=
                    entry.palinstrophy_normalization *
                    matrix.power_one_scale;
            }
        }
        if (matrix.selected_enstrophy > 0.0L &&
            matrix.selected_palinstrophy > 0.0L) {
            const SpectralReal cauchy_factor =
                1.5L * std::abs(matrix.power_one_scale) /
                matrix.selected_palinstrophy;
            const SpectralReal z = matrix.selected_enstrophy;
            const SpectralReal p = matrix.selected_palinstrophy;
            const SpectralReal core_h1 = std::max(
                row.core_aggregate_h1_norm2, 0.0L);
            const SpectralReal core_h2 = std::max(
                row.core_aggregate_h2_norm2, 0.0L);
            const SpectralReal tail_h1 = std::max(
                row.tail_aggregate_h1_norm2, 0.0L);
            const SpectralReal tail_h2 = std::max(
                row.tail_aggregate_h2_norm2, 0.0L);
            row.core_stretching_tail_cross_cauchy_bound =
                cauchy_factor * std::sqrt(z * core_h1 * p * tail_h2);
            row.tail_stretching_core_cross_cauchy_bound =
                cauchy_factor * std::sqrt(z * tail_h1 * p * core_h2);
            row.tail_stretching_tail_cross_cauchy_bound =
                cauchy_factor * std::sqrt(z * tail_h1 * p * tail_h2);
            auto alignment = [](SpectralReal value,
                                SpectralReal norm_product) {
                return norm_product > 0.0L
                    ? std::abs(value) / std::sqrt(norm_product)
                    : 0.0L;
            };
            row.core_stretching_h1_alignment = alignment(
                row.core_stretching, z * core_h1);
            row.tail_stretching_h1_alignment = alignment(
                row.tail_stretching, z * tail_h1);
            row.core_palinstrophy_cross_h2_alignment = alignment(
                row.core_palinstrophy_cross, p * core_h2);
            row.tail_palinstrophy_cross_h2_alignment = alignment(
                row.tail_palinstrophy_cross, p * tail_h2);
        }
        auto cauchy_ratio = [](SpectralReal actual, SpectralReal bound) {
            return bound > 0.0L ? std::abs(actual) / bound : 0.0L;
        };
        row.core_stretching_tail_cross_cauchy_ratio = cauchy_ratio(
            row.core_stretching_tail_cross_power_one,
            row.core_stretching_tail_cross_cauchy_bound);
        row.tail_stretching_core_cross_cauchy_ratio = cauchy_ratio(
            row.tail_stretching_core_cross_power_one,
            row.tail_stretching_core_cross_cauchy_bound);
        row.tail_stretching_tail_cross_cauchy_ratio = cauchy_ratio(
            row.tail_stretching_tail_cross_power_one,
            row.tail_stretching_tail_cross_cauchy_bound);
        row.joint_cross_tail_cauchy_bound =
            row.core_stretching_tail_cross_cauchy_bound +
            row.tail_stretching_core_cross_cauchy_bound +
            row.tail_stretching_tail_cross_cauchy_bound;
        row.joint_cross_tail_cauchy_ratio = cauchy_ratio(
            row.open_palinstrophy_normalization_power_one,
            row.joint_cross_tail_cauchy_bound);
        row.maximum_alignment_product_reconstruction_error = std::max({
            relative_error(
                row.core_stretching_tail_cross_cauchy_ratio,
                row.core_stretching_h1_alignment *
                    row.tail_palinstrophy_cross_h2_alignment),
            relative_error(
                row.tail_stretching_core_cross_cauchy_ratio,
                row.tail_stretching_h1_alignment *
                    row.core_palinstrophy_cross_h2_alignment),
            relative_error(
                row.tail_stretching_tail_cross_cauchy_ratio,
                row.tail_stretching_h1_alignment *
                    row.tail_palinstrophy_cross_h2_alignment)});
        row.open_power_one = row.core_tail_power_one +
            row.tail_internal_power_one;
        row.reconstructed_power_one = row.core_internal_power_one +
            row.open_power_one;
        if (open_square_sum > 0.0L) {
            row.open_effective_height_pairs =
                row.open_absolute_power_one_sum *
                row.open_absolute_power_one_sum / open_square_sum;
        }
        if (row.open_absolute_power_one_sum > 0.0L) {
            row.open_signed_alignment =
                std::abs(row.open_power_one) /
                row.open_absolute_power_one_sum;
            for (const auto& entry : matrix.entries) {
                if (entry.second_shell > row.last_core_shell) {
                    row.dominant_open_pair_fraction = std::max(
                        row.dominant_open_pair_fraction,
                        std::abs(entry.power_one) /
                            row.open_absolute_power_one_sum);
                }
            }
        }
        row.reconstruction_error = relative_error(
            row.reconstructed_power_one, matrix.selected_power_one);
        row.component_reconstruction_error = relative_error(
            row.open_outer_square_power_one +
                row.open_advected_commutator_power_one +
                row.open_advecting_nested_power_one +
                row.open_enstrophy_normalization_power_one +
                row.open_palinstrophy_normalization_power_one,
            row.open_power_one);
        row.palinstrophy_factorization_error = relative_error(
            row.core_stretching_tail_cross_power_one +
                row.tail_stretching_core_cross_power_one +
                row.tail_stretching_tail_cross_power_one,
            row.open_palinstrophy_normalization_power_one);
        report.maximum_reconstruction_error = std::max(
            report.maximum_reconstruction_error,
            row.reconstruction_error);
        report.maximum_component_reconstruction_error = std::max(
            report.maximum_component_reconstruction_error,
            row.component_reconstruction_error);
        report.maximum_palinstrophy_factorization_error = std::max(
            report.maximum_palinstrophy_factorization_error,
            row.palinstrophy_factorization_error);
        report.maximum_alignment_product_reconstruction_error = std::max(
            report.maximum_alignment_product_reconstruction_error,
            row.maximum_alignment_product_reconstruction_error);
        report.maximum_normalization_cauchy_ratio = std::max({
            report.maximum_normalization_cauchy_ratio,
            row.core_stretching_tail_cross_cauchy_ratio,
            row.tail_stretching_core_cross_cauchy_ratio,
            row.tail_stretching_tail_cross_cauchy_ratio,
            row.joint_cross_tail_cauchy_ratio});
        report.rows.push_back(row);
    }
    report.exact_cumulative_decomposition =
        report.maximum_reconstruction_error < 1e-13L &&
        report.maximum_component_reconstruction_error < 1e-13L &&
        report.maximum_palinstrophy_factorization_error < 1e-13L &&
        report.maximum_alignment_product_reconstruction_error < 1e-13L;
    report.finite_normalization_cauchy_inequalities_verified =
        report.maximum_normalization_cauchy_ratio <= 1.0L + 1e-13L;
    return report;
}

}  // namespace lemma
