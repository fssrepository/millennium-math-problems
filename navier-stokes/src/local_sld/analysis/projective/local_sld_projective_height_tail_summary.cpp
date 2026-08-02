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
        SpectralReal open_square_sum = 0.0L;
        for (const auto& entry : matrix.entries) {
            const bool first_core =
                entry.first_shell <= row.last_core_shell;
            const bool second_core =
                entry.second_shell <= row.last_core_shell;
            if (first_core && second_core) {
                row.core_internal_power_one += entry.power_one;
            } else if (first_core) {
                row.core_tail_power_one += entry.power_one;
                row.open_absolute_power_one_sum +=
                    std::abs(entry.power_one);
                open_square_sum += entry.power_one * entry.power_one;
            } else {
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
        report.maximum_reconstruction_error = std::max(
            report.maximum_reconstruction_error,
            row.reconstruction_error);
        report.maximum_component_reconstruction_error = std::max(
            report.maximum_component_reconstruction_error,
            row.component_reconstruction_error);
        report.rows.push_back(row);
    }
    report.exact_cumulative_decomposition =
        report.maximum_reconstruction_error < 1e-13L &&
        report.maximum_component_reconstruction_error < 1e-13L;
    return report;
}

}  // namespace lemma
