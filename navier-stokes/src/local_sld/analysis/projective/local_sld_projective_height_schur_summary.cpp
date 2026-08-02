#include "local_sld_projective_height_schur_summary.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {

LocalSldProjectiveHeightSchurReport
LocalSldProjectiveHeightSchurSummary::summarize(
    const LocalSldProjectiveHeightMatrixReport& matrix) {
    if (matrix.shells.empty()) {
        throw std::invalid_argument(
            "projective height Schur summary requires height shells");
    }
    const std::size_t shell_count = matrix.shells.size();
    std::vector<SpectralReal> diagonal_envelope(shell_count, 0.0L);
    std::vector<SpectralReal> commutator_paired_diagonal(
        shell_count, 0.0L);
    std::vector<SpectralReal> outer_power_diagonal(
        shell_count, 0.0L);
    std::vector<SpectralReal> dynamic_response_diagonal(
        shell_count, 0.0L);
    for (const auto& entry : matrix.entries) {
        if (entry.first_shell == entry.second_shell) {
            diagonal_envelope[
                static_cast<std::size_t>(entry.first_shell)] =
                entry.absolute_component_power_one_envelope;
            commutator_paired_diagonal[
                static_cast<std::size_t>(entry.first_shell)] =
                entry.commutator_paired_power_one_envelope;
            outer_power_diagonal[
                static_cast<std::size_t>(entry.first_shell)] =
                std::abs(entry.outer_square * matrix.power_one_scale);
            const auto& shell = matrix.shells[
                static_cast<std::size_t>(entry.first_shell)];
            dynamic_response_diagonal[
                static_cast<std::size_t>(entry.first_shell)] =
                std::abs(matrix.power_one_scale) *
                (shell.aggregate_h1_norm2 +
                 shell.dynamic_response_hminus1_norm2);
        }
    }
    LocalSldProjectiveHeightSchurReport report;
    report.gaps.resize(shell_count);
    for (std::size_t gap = 0; gap < shell_count; ++gap) {
        report.gaps[gap].shell_gap = static_cast<int>(gap);
    }
    std::vector<SpectralReal> row_sums(shell_count, 0.0L);
    std::vector<SpectralReal> commutator_paired_row_sums(
        shell_count, 0.0L);
    std::vector<SpectralReal> commutator_paired_outer_row_sums(
        shell_count, 0.0L);
    std::vector<SpectralReal> dynamic_paired_outer_row_sums(
        shell_count, 0.0L);
    std::vector<SpectralReal> dynamic_paired_response_row_sums(
        shell_count, 0.0L);
    for (const SpectralReal weight : outer_power_diagonal) {
        report.commutator_paired_outer_weight += weight;
    }
    for (const SpectralReal weight : dynamic_response_diagonal) {
        report.dynamic_paired_response_weight += weight;
    }
    for (const auto& source : matrix.entries) {
        LocalSldProjectiveHeightSchurEntry entry;
        entry.first_shell = source.first_shell;
        entry.second_shell = source.second_shell;
        entry.shell_gap = source.second_shell - source.first_shell;
        entry.power_one = source.power_one;
        entry.absolute_component_power_one_envelope =
            source.absolute_component_power_one_envelope;
        report.total_component_envelope +=
            entry.absolute_component_power_one_envelope;
        if (entry.first_shell == entry.second_shell) {
            report.diagonal_component_envelope +=
                entry.absolute_component_power_one_envelope;
        }
        report.commutator_paired_total_envelope +=
            source.commutator_paired_power_one_envelope;
        report.dynamic_paired_total_envelope +=
            source.dynamic_paired_power_one_envelope;
        if (entry.first_shell == entry.second_shell) {
            report.commutator_paired_diagonal_envelope +=
                source.commutator_paired_power_one_envelope;
        }
        const SpectralReal paired_scale = std::sqrt(
            commutator_paired_diagonal[
                static_cast<std::size_t>(entry.first_shell)] *
            commutator_paired_diagonal[
                static_cast<std::size_t>(entry.second_shell)]);
        if (paired_scale > 1e-30L) {
            const SpectralReal symmetry_factor =
                entry.first_shell == entry.second_shell
                ? 1.0L : 2.0L;
            const SpectralReal ratio =
                source.commutator_paired_power_one_envelope /
                (symmetry_factor * paired_scale);
            commutator_paired_row_sums[
                static_cast<std::size_t>(entry.first_shell)] += ratio;
            if (entry.first_shell != entry.second_shell) {
                commutator_paired_row_sums[
                    static_cast<std::size_t>(entry.second_shell)] += ratio;
            }
        } else if (entry.first_shell != entry.second_shell &&
                   source.commutator_paired_power_one_envelope > 1e-30L) {
            ++report
                .commutator_paired_unscaled_off_diagonal_pair_count;
        }
        const SpectralReal paired_outer_scale = std::sqrt(
            outer_power_diagonal[
                static_cast<std::size_t>(entry.first_shell)] *
            outer_power_diagonal[
                static_cast<std::size_t>(entry.second_shell)]);
        SpectralReal paired_outer_ratio = 0.0L;
        SpectralReal dynamic_outer_ratio = 0.0L;
        const bool has_paired_outer_scale =
            paired_outer_scale > 1e-30L;
        if (has_paired_outer_scale) {
            const SpectralReal symmetry_factor =
                entry.first_shell == entry.second_shell
                ? 1.0L : 2.0L;
            paired_outer_ratio =
                source.commutator_paired_power_one_envelope /
                (symmetry_factor * paired_outer_scale);
            commutator_paired_outer_row_sums[
                static_cast<std::size_t>(entry.first_shell)] +=
                paired_outer_ratio;
            if (entry.first_shell != entry.second_shell) {
                commutator_paired_outer_row_sums[
                    static_cast<std::size_t>(entry.second_shell)] +=
                    paired_outer_ratio;
            }
            dynamic_outer_ratio =
                source.dynamic_paired_power_one_envelope /
                (symmetry_factor * paired_outer_scale);
            dynamic_paired_outer_row_sums[
                static_cast<std::size_t>(entry.first_shell)] +=
                dynamic_outer_ratio;
            if (entry.first_shell != entry.second_shell) {
                dynamic_paired_outer_row_sums[
                    static_cast<std::size_t>(entry.second_shell)] +=
                    dynamic_outer_ratio;
            }
        } else if (entry.first_shell != entry.second_shell &&
                   source.commutator_paired_power_one_envelope > 1e-30L) {
            ++report
                .commutator_paired_outer_unscaled_off_diagonal_pair_count;
        }
        if (!has_paired_outer_scale &&
            entry.first_shell != entry.second_shell &&
            source.dynamic_paired_power_one_envelope > 1e-30L) {
            ++report.dynamic_paired_outer_unscaled_off_diagonal_pair_count;
        }
        const SpectralReal response_scale = std::sqrt(
            dynamic_response_diagonal[
                static_cast<std::size_t>(entry.first_shell)] *
            dynamic_response_diagonal[
                static_cast<std::size_t>(entry.second_shell)]);
        SpectralReal dynamic_response_ratio = 0.0L;
        const bool has_response_scale = response_scale > 1e-30L;
        if (has_response_scale) {
            const SpectralReal symmetry_factor =
                entry.first_shell == entry.second_shell
                ? 1.0L : 2.0L;
            dynamic_response_ratio =
                source.dynamic_paired_power_one_envelope /
                (symmetry_factor * response_scale);
            dynamic_paired_response_row_sums[
                static_cast<std::size_t>(entry.first_shell)] +=
                dynamic_response_ratio;
            if (entry.first_shell != entry.second_shell) {
                dynamic_paired_response_row_sums[
                    static_cast<std::size_t>(entry.second_shell)] +=
                    dynamic_response_ratio;
            }
        } else if (entry.first_shell != entry.second_shell &&
                   source.dynamic_paired_power_one_envelope > 1e-30L) {
            ++report
                .dynamic_paired_response_unscaled_off_diagonal_pair_count;
        }
        const SpectralReal first_diagonal = diagonal_envelope[
            static_cast<std::size_t>(entry.first_shell)];
        const SpectralReal second_diagonal = diagonal_envelope[
            static_cast<std::size_t>(entry.second_shell)];
        entry.geometric_diagonal_scale = std::sqrt(
            first_diagonal * second_diagonal);
        entry.has_nonzero_diagonal_scale =
            entry.geometric_diagonal_scale > 1e-30L;
        if (entry.has_nonzero_diagonal_scale) {
            const SpectralReal symmetry_factor =
                entry.first_shell == entry.second_shell
                ? 1.0L : 2.0L;
            entry.symmetric_geometric_ratio =
                entry.absolute_component_power_one_envelope /
                (symmetry_factor * entry.geometric_diagonal_scale);
            report.maximum_symmetric_geometric_ratio = std::max(
                report.maximum_symmetric_geometric_ratio,
                entry.symmetric_geometric_ratio);
            row_sums[static_cast<std::size_t>(entry.first_shell)] +=
                entry.symmetric_geometric_ratio;
            if (entry.first_shell != entry.second_shell) {
                row_sums[static_cast<std::size_t>(entry.second_shell)] +=
                    entry.symmetric_geometric_ratio;
            }
        } else if (entry.first_shell != entry.second_shell &&
                   entry.absolute_component_power_one_envelope > 1e-30L) {
            ++report.unscaled_off_diagonal_pair_count;
        }
        auto& gap = report.gaps[static_cast<std::size_t>(entry.shell_gap)];
        ++gap.pair_count;
        gap.signed_power_one_sum += entry.power_one;
        gap.absolute_power_one_sum += std::abs(entry.power_one);
        gap.absolute_component_envelope_sum +=
            entry.absolute_component_power_one_envelope;
        gap.commutator_paired_envelope_sum +=
            source.commutator_paired_power_one_envelope;
        gap.dynamic_paired_envelope_sum +=
            source.dynamic_paired_power_one_envelope;
        const SpectralReal absolute_scale =
            std::abs(matrix.power_one_scale);
        const SpectralReal commutator_envelope =
            std::abs(source.outer_square +
                     source.advected_commutator) * absolute_scale;
        const SpectralReal remainder_envelope =
            (std::abs(source.advecting_nested) +
             std::abs(source.enstrophy_normalization) +
             std::abs(source.palinstrophy_normalization)) *
            absolute_scale;
        gap.commutator_term_envelope_sum += commutator_envelope;
        gap.remainder_terms_envelope_sum += remainder_envelope;
        if (has_paired_outer_scale) {
            gap.commutator_paired_outer_maximum_symmetric_geometric_ratio =
                std::max(
                    gap
                        .commutator_paired_outer_maximum_symmetric_geometric_ratio,
                    paired_outer_ratio);
            const SpectralReal symmetry_factor =
                entry.first_shell == entry.second_shell
                ? 1.0L : 2.0L;
            gap.commutator_outer_maximum_symmetric_geometric_ratio =
                std::max(
                    gap
                        .commutator_outer_maximum_symmetric_geometric_ratio,
                    commutator_envelope /
                        (symmetry_factor * paired_outer_scale));
            gap.remainder_outer_maximum_symmetric_geometric_ratio =
                std::max(
                    gap.remainder_outer_maximum_symmetric_geometric_ratio,
                    remainder_envelope /
                        (symmetry_factor * paired_outer_scale));
            gap.dynamic_paired_outer_maximum_symmetric_geometric_ratio =
                std::max(
                    gap
                        .dynamic_paired_outer_maximum_symmetric_geometric_ratio,
                    dynamic_outer_ratio);
        } else if (entry.first_shell != entry.second_shell &&
                   source.commutator_paired_power_one_envelope > 1e-30L) {
            ++gap.commutator_paired_outer_unscaled_pair_count;
        }
        if (!has_paired_outer_scale &&
            entry.first_shell != entry.second_shell &&
            source.dynamic_paired_power_one_envelope > 1e-30L) {
            ++gap.dynamic_paired_outer_unscaled_pair_count;
        }
        if (has_response_scale) {
            gap
                .dynamic_paired_response_maximum_symmetric_geometric_ratio =
                std::max(
                    gap
                        .dynamic_paired_response_maximum_symmetric_geometric_ratio,
                    dynamic_response_ratio);
        } else if (entry.first_shell != entry.second_shell &&
                   source.dynamic_paired_power_one_envelope > 1e-30L) {
            ++gap.dynamic_paired_response_unscaled_pair_count;
        }
        if (entry.has_nonzero_diagonal_scale) {
            gap.maximum_symmetric_geometric_ratio = std::max(
                gap.maximum_symmetric_geometric_ratio,
                entry.symmetric_geometric_ratio);
        } else if (entry.first_shell != entry.second_shell &&
                   entry.absolute_component_power_one_envelope > 1e-30L) {
            ++gap.unscaled_pair_count;
        }
        report.entries.push_back(entry);
    }
    report.maximum_weighted_row_sum = *std::max_element(
        row_sums.begin(), row_sums.end());
    report.commutator_paired_maximum_weighted_row_sum =
        *std::max_element(
            commutator_paired_row_sums.begin(),
            commutator_paired_row_sums.end());
    report.commutator_paired_outer_maximum_weighted_row_sum =
        *std::max_element(
            commutator_paired_outer_row_sums.begin(),
            commutator_paired_outer_row_sums.end());
    report.dynamic_paired_outer_maximum_weighted_row_sum =
        *std::max_element(
            dynamic_paired_outer_row_sums.begin(),
            dynamic_paired_outer_row_sums.end());
    report.dynamic_paired_response_maximum_weighted_row_sum =
        *std::max_element(
            dynamic_paired_response_row_sums.begin(),
            dynamic_paired_response_row_sums.end());
    for (const auto& gap : report.gaps) {
        report.dynamic_paired_outer_gap_ratio_sum +=
            gap.dynamic_paired_outer_maximum_symmetric_geometric_ratio;
        report.dynamic_paired_outer_gap_one_decay_constant = std::max(
            report.dynamic_paired_outer_gap_one_decay_constant,
            std::ldexp(
                gap.dynamic_paired_outer_maximum_symmetric_geometric_ratio,
                gap.shell_gap));
        report.dynamic_paired_response_gap_ratio_sum +=
            gap
                .dynamic_paired_response_maximum_symmetric_geometric_ratio;
        report.dynamic_paired_response_gap_one_decay_constant = std::max(
            report.dynamic_paired_response_gap_one_decay_constant,
            std::ldexp(
                gap
                    .dynamic_paired_response_maximum_symmetric_geometric_ratio,
                gap.shell_gap));
    }
    report.weighted_schur_upper_bound =
        report.maximum_weighted_row_sum *
        report.diagonal_component_envelope;
    if (report.weighted_schur_upper_bound > 0.0L) {
        report.upper_bound_ratio =
            report.total_component_envelope /
            report.weighted_schur_upper_bound;
    }
    report.commutator_paired_weighted_schur_upper_bound =
        report.commutator_paired_maximum_weighted_row_sum *
        report.commutator_paired_diagonal_envelope;
    if (report.commutator_paired_weighted_schur_upper_bound > 0.0L) {
        report.commutator_paired_upper_bound_ratio =
            report.commutator_paired_total_envelope /
            report.commutator_paired_weighted_schur_upper_bound;
    }
    report.commutator_paired_outer_weighted_schur_upper_bound =
        report.commutator_paired_outer_maximum_weighted_row_sum *
        report.commutator_paired_outer_weight;
    if (report.commutator_paired_outer_weighted_schur_upper_bound >
        0.0L) {
        report.commutator_paired_outer_upper_bound_ratio =
            report.commutator_paired_total_envelope /
            report.commutator_paired_outer_weighted_schur_upper_bound;
    }
    report.dynamic_paired_outer_weighted_schur_upper_bound =
        report.dynamic_paired_outer_maximum_weighted_row_sum *
        report.commutator_paired_outer_weight;
    if (report.dynamic_paired_outer_weighted_schur_upper_bound > 0.0L) {
        report.dynamic_paired_outer_upper_bound_ratio =
            report.dynamic_paired_total_envelope /
            report.dynamic_paired_outer_weighted_schur_upper_bound;
    }
    report.dynamic_paired_response_weighted_schur_upper_bound =
        report.dynamic_paired_response_maximum_weighted_row_sum *
        report.dynamic_paired_response_weight;
    if (report.dynamic_paired_response_weighted_schur_upper_bound > 0.0L) {
        report.dynamic_paired_response_upper_bound_ratio =
            report.dynamic_paired_total_envelope /
            report.dynamic_paired_response_weighted_schur_upper_bound;
    }
    std::vector<SpectralReal> outer_diagonal(shell_count, 0.0L);
    for (const auto& entry : matrix.entries) {
        if (entry.first_shell == entry.second_shell) {
            outer_diagonal[
                static_cast<std::size_t>(entry.first_shell)] =
                std::abs(entry.outer_square);
        }
    }
    auto analyze_component = [&](const char* name, auto value) {
        LocalSldProjectiveHeightComponentSchurRow component;
        component.component = name;
        std::vector<SpectralReal> component_diagonal(
            shell_count, 0.0L);
        for (const auto& entry : matrix.entries) {
            const SpectralReal scaled_magnitude =
                std::abs(value(entry) * matrix.power_one_scale);
            component.total_power_one_envelope += scaled_magnitude;
            if (entry.first_shell == entry.second_shell) {
                component.diagonal_power_one_envelope +=
                    scaled_magnitude;
                component_diagonal[
                    static_cast<std::size_t>(entry.first_shell)] =
                    std::abs(value(entry));
            }
        }
        std::vector<SpectralReal> component_row_sums(
            shell_count, 0.0L);
        for (const auto& entry : matrix.entries) {
            const SpectralReal first = component_diagonal[
                static_cast<std::size_t>(entry.first_shell)];
            const SpectralReal second = component_diagonal[
                static_cast<std::size_t>(entry.second_shell)];
            const SpectralReal scale = std::sqrt(first * second);
            const SpectralReal magnitude = std::abs(value(entry));
            if (scale > 1e-30L) {
                const SpectralReal symmetry_factor =
                    entry.first_shell == entry.second_shell
                    ? 1.0L : 2.0L;
                const SpectralReal ratio =
                    magnitude / (symmetry_factor * scale);
                component_row_sums[
                    static_cast<std::size_t>(entry.first_shell)] += ratio;
                if (entry.first_shell != entry.second_shell) {
                    component_row_sums[
                        static_cast<std::size_t>(entry.second_shell)] +=
                        ratio;
                    component.maximum_off_diagonal_geometric_ratio =
                        std::max(
                            component
                                .maximum_off_diagonal_geometric_ratio,
                            ratio);
                }
            } else if (entry.first_shell != entry.second_shell &&
                       magnitude > 1e-30L) {
                ++component.unscaled_off_diagonal_pair_count;
            }
        }
        component.maximum_weighted_row_sum = *std::max_element(
            component_row_sums.begin(), component_row_sums.end());
        std::vector<SpectralReal> outer_row_sums(
            shell_count, 0.0L);
        for (const auto& entry : matrix.entries) {
            const SpectralReal first = outer_diagonal[
                static_cast<std::size_t>(entry.first_shell)];
            const SpectralReal second = outer_diagonal[
                static_cast<std::size_t>(entry.second_shell)];
            const SpectralReal scale = std::sqrt(first * second);
            const SpectralReal magnitude = std::abs(value(entry));
            if (scale > 1e-30L) {
                const SpectralReal symmetry_factor =
                    entry.first_shell == entry.second_shell
                    ? 1.0L : 2.0L;
                const SpectralReal ratio =
                    magnitude / (symmetry_factor * scale);
                outer_row_sums[
                    static_cast<std::size_t>(entry.first_shell)] += ratio;
                if (entry.first_shell == entry.second_shell) {
                    component.maximum_outer_weighted_diagonal_ratio =
                        std::max(
                            component
                                .maximum_outer_weighted_diagonal_ratio,
                            ratio);
                } else {
                    outer_row_sums[
                        static_cast<std::size_t>(entry.second_shell)] +=
                        ratio;
                    component
                        .maximum_outer_weighted_off_diagonal_ratio =
                        std::max(
                            component
                                .maximum_outer_weighted_off_diagonal_ratio,
                            ratio);
                }
            } else if (magnitude > 1e-30L) {
                ++component.outer_unscaled_pair_count;
            }
        }
        component.maximum_outer_weighted_row_sum = *std::max_element(
            outer_row_sums.begin(), outer_row_sums.end());
        report.components.push_back(component);
    };
    analyze_component("outer_square", [](const auto& entry) {
        return entry.outer_square;
    });
    analyze_component("advected_commutator", [](const auto& entry) {
        return entry.advected_commutator;
    });
    analyze_component("advecting_nested", [](const auto& entry) {
        return entry.advecting_nested;
    });
    analyze_component("enstrophy_normalization", [](const auto& entry) {
        return entry.enstrophy_normalization;
    });
    analyze_component("palinstrophy_normalization", [](const auto& entry) {
        return entry.palinstrophy_normalization;
    });
    report.finite_matrix_exact =
        matrix.exact_height_matrix_decomposition;
    report.finite_schur_inequality_verified =
        report.unscaled_off_diagonal_pair_count == 0 &&
        report.total_component_envelope <=
            report.weighted_schur_upper_bound * (1.0L + 1e-14L);
    report.finite_commutator_paired_schur_inequality_verified =
        report.commutator_paired_unscaled_off_diagonal_pair_count == 0 &&
        report.commutator_paired_total_envelope <=
            report.commutator_paired_weighted_schur_upper_bound *
                (1.0L + 1e-14L);
    report.finite_commutator_paired_outer_schur_inequality_verified =
        report
            .commutator_paired_outer_unscaled_off_diagonal_pair_count == 0 &&
        report.commutator_paired_total_envelope <=
            report.commutator_paired_outer_weighted_schur_upper_bound *
                (1.0L + 1e-14L);
    report.finite_dynamic_paired_outer_schur_inequality_verified =
        report.dynamic_paired_outer_unscaled_off_diagonal_pair_count == 0 &&
        report.dynamic_paired_total_envelope <=
            report.dynamic_paired_outer_weighted_schur_upper_bound *
                (1.0L + 1e-14L);
    report.finite_dynamic_paired_response_schur_inequality_verified =
        report
            .dynamic_paired_response_unscaled_off_diagonal_pair_count == 0 &&
        report.dynamic_paired_total_envelope <=
            report.dynamic_paired_response_weighted_schur_upper_bound *
                (1.0L + 1e-14L);
    return report;
}

}  // namespace lemma
