#include "local_sld_projective_normalization_tail_schur_ledger.hpp"

#include "local_sld_projective_normalization_cauchy_objective.hpp"
#include "projective_advection_decomposition.hpp"
#include "projective_height_shell_partition.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace lemma {
namespace {

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "normalization tail Schur pairing layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "normalization tail Schur Laplacian layout mismatch");
    }
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

void add(
    SpectralIncrement& target,
    const SpectralIncrement& source) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "normalization tail Schur sum layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] += source[mode][coordinate];
        }
    }
}

SpectralReal relative_error(
    SpectralReal computed,
    SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

struct ShellEvaluation {
    LocalSldProjectiveNormalizationTailSchurShell row;
    std::vector<std::size_t> group_indices;
    SpectralIncrement ab;
};

}  // namespace

LocalSldProjectiveNormalizationTailSchurReport
LocalSldProjectiveNormalizationTailSchurLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads) {
    if (core_maximum_height < 1) {
        throw std::invalid_argument(
            "normalization tail Schur core height must be positive");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "normalization tail Schur threads must be 1..256");
    }
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    LocalSldProjectiveNormalizationTailSchurReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.core_maximum_height = core_maximum_height;
    report.selected_shape_count = groups.size();

    std::map<int, ShellEvaluation> tail_shells;
    for (std::size_t index = 0; index < groups.size(); ++index) {
        const SpectralInteger height =
            groups[index].primitive_squared_lengths[2];
        if (height <= core_maximum_height) {
            continue;
        }
        ++report.tail_shape_count;
        const int shell = ProjectiveHeightShellPartition::shell_index(
            height);
        ShellEvaluation& target = tail_shells[shell];
        target.row.shell = shell;
        target.row.minimum_primitive_height =
            target.row.minimum_primitive_height == 0
            ? height
            : std::min(target.row.minimum_primitive_height, height);
        target.row.maximum_primitive_height = std::max(
            target.row.maximum_primitive_height, height);
        ++target.row.shape_count;
        target.row.interaction_count += groups[index].interactions.size();
        target.group_indices.push_back(index);
    }

    std::vector<ShellEvaluation> shells;
    shells.reserve(tail_shells.size());
    for (auto& [shell, evaluation] : tail_shells) {
        static_cast<void>(shell);
        SpectralIncrement b =
            ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
                state, groups, evaluation.group_indices,
                state.velocity, state.velocity, threads);
        evaluation.ab = laplacian_weight(state, b);
        evaluation.row.aggregate_h2_norm2 = pairing(
            evaluation.ab, evaluation.ab);
        report.diagonal_tail_h2_norm2 +=
            evaluation.row.aggregate_h2_norm2;
        shells.push_back(std::move(evaluation));
    }

    SpectralIncrement aggregate_ab(state.waves.size());
    for (const ShellEvaluation& shell : shells) {
        add(aggregate_ab, shell.ab);
    }
    report.aggregate_tail_h2_norm2 = pairing(
        aggregate_ab, aggregate_ab);

    std::vector<SpectralReal> normalized_rows(shells.size(), 0.0L);
    std::map<int, LocalSldProjectiveNormalizationTailSchurGap> gaps;
    for (std::size_t first = 0; first < shells.size(); ++first) {
        for (std::size_t second = first; second < shells.size(); ++second) {
            const SpectralReal gram = pairing(
                shells[first].ab, shells[second].ab);
            const SpectralReal multiplicity =
                first == second ? 1.0L : 2.0L;
            report.reconstructed_tail_h2_norm2 += multiplicity * gram;
            report.absolute_gram_tail_h2 +=
                multiplicity * std::abs(gram);
            const SpectralReal denominator = std::sqrt(std::max(
                shells[first].row.aggregate_h2_norm2 *
                    shells[second].row.aggregate_h2_norm2,
                0.0L));
            if (denominator > 0.0L) {
                const SpectralReal normalized =
                    std::abs(gram) / denominator;
                const int gap = shells[second].row.shell -
                    shells[first].row.shell;
                const SpectralReal half_decay_weighted = std::ldexp(
                    normalized, gap);
                normalized_rows[first] += normalized;
                if (second != first) {
                    normalized_rows[second] += normalized;
                }
                report.pairs.push_back({
                    shells[first].row.shell,
                    shells[second].row.shell,
                    gap,
                    gram,
                    normalized,
                    half_decay_weighted});
                auto& gap_row = gaps[gap];
                gap_row.shell_gap = gap;
                ++gap_row.pair_count;
                gap_row.maximum_normalized_absolute_gram = std::max(
                    gap_row.maximum_normalized_absolute_gram,
                    normalized);
                gap_row.sum_normalized_absolute_gram += normalized;
                gap_row.maximum_half_decay_weighted_correlation = std::max(
                    gap_row.maximum_half_decay_weighted_correlation,
                    half_decay_weighted);
                if (gap > 0) {
                    report.maximum_half_decay_weighted_correlation =
                        std::max(
                            report.maximum_half_decay_weighted_correlation,
                            half_decay_weighted);
                }
            }
        }
    }
    for (auto& [gap, row] : gaps) {
        static_cast<void>(gap);
        report.gaps.push_back(std::move(row));
    }
    for (std::size_t index = 0; index < shells.size(); ++index) {
        shells[index].row.normalized_absolute_gram_row_sum =
            normalized_rows[index];
        if (report.diagonal_tail_h2_norm2 > 0.0L) {
            shells[index].row.diagonal_tail_fraction =
                shells[index].row.aggregate_h2_norm2 /
                report.diagonal_tail_h2_norm2;
        }
        report.maximum_normalized_absolute_gram_row_sum = std::max(
            report.maximum_normalized_absolute_gram_row_sum,
            normalized_rows[index]);
        report.shells.push_back(std::move(shells[index].row));
    }
    report.schur_tail_h2_upper_bound =
        report.maximum_normalized_absolute_gram_row_sum *
        report.diagonal_tail_h2_norm2;
    report.half_decay_implied_row_bound =
        1.0L + 2.0L *
            report.maximum_half_decay_weighted_correlation;
    report.half_decay_tail_h2_upper_bound =
        report.half_decay_implied_row_bound *
        report.diagonal_tail_h2_norm2;
    if (report.diagonal_tail_h2_norm2 > 0.0L) {
        report.aggregate_to_diagonal_ratio =
            report.aggregate_tail_h2_norm2 /
            report.diagonal_tail_h2_norm2;
    }
    if (report.absolute_gram_tail_h2 > 0.0L) {
        report.aggregate_to_absolute_gram_ratio =
            report.aggregate_tail_h2_norm2 /
            report.absolute_gram_tail_h2;
    }
    if (report.schur_tail_h2_upper_bound > 0.0L) {
        report.aggregate_to_schur_bound_ratio =
            report.aggregate_tail_h2_norm2 /
            report.schur_tail_h2_upper_bound;
    }
    report.tail_h2_reconstruction_error = relative_error(
        report.reconstructed_tail_h2_norm2,
        report.aggregate_tail_h2_norm2);

    const auto cauchy = LocalSldProjectiveNormalizationCauchyObjective(
        dynamics, selection, core_maximum_height, threads).evaluate(state);
    report.full_stretching = cauchy.full_stretching;
    report.enstrophy = cauchy.enstrophy;
    report.palinstrophy = cauchy.palinstrophy;
    report.selected_aggregate_h1_norm2 =
        cauchy.selected_aggregate_h1_norm2;
    report.objective_tail_h2_reconstruction_error = relative_error(
        report.aggregate_tail_h2_norm2,
        cauchy.tail_aggregate_h2_norm2);
    const SpectralReal denominator =
        cauchy.enstrophy * cauchy.enstrophy * cauchy.enstrophy *
        cauchy.palinstrophy * cauchy.palinstrophy *
        cauchy.palinstrophy * cauchy.palinstrophy * cauchy.palinstrophy;
    SpectralReal common_factor = 0.0L;
    if (denominator > 1e-60L) {
        common_factor = 2.25L * cauchy.full_stretching *
            cauchy.full_stretching *
            cauchy.selected_aggregate_h1_norm2 / denominator;
    }
    report.normalization_common_factor = common_factor;
    report.squared_cauchy_majorant =
        common_factor * report.aggregate_tail_h2_norm2;
    report.diagonal_squared_majorant =
        common_factor * report.diagonal_tail_h2_norm2;
    report.schur_squared_majorant =
        common_factor * report.schur_tail_h2_upper_bound;
    const SpectralReal half_decay_squared_majorant =
        common_factor * report.half_decay_tail_h2_upper_bound;
    const SpectralReal height_factor = std::sqrt(
        static_cast<SpectralReal>(core_maximum_height));
    report.height_half_compensated_squared_cauchy_majorant =
        height_factor * report.squared_cauchy_majorant;
    report.height_half_compensated_diagonal_squared_majorant =
        height_factor * report.diagonal_squared_majorant;
    report.height_half_compensated_schur_squared_majorant =
        height_factor * report.schur_squared_majorant;
    report.height_half_compensated_half_decay_squared_majorant =
        height_factor * half_decay_squared_majorant;

    const SpectralReal tolerance = 2e-12L;
    report.exact_tail_h2_reconstruction =
        report.tail_h2_reconstruction_error < tolerance &&
        report.objective_tail_h2_reconstruction_error < tolerance;
    report.finite_schur_inequality_verified =
        report.aggregate_tail_h2_norm2 <=
            (1.0L + tolerance) * report.schur_tail_h2_upper_bound;
    report.finite_half_gap_decay_inequality_verified =
        report.maximum_normalized_absolute_gram_row_sum <=
            (1.0L + tolerance) * report.half_decay_implied_row_bound &&
        report.aggregate_tail_h2_norm2 <=
            (1.0L + tolerance) *
                report.half_decay_tail_h2_upper_bound;
    report.finite = cauchy.finite &&
        std::isfinite(report.schur_squared_majorant) &&
        report.exact_tail_h2_reconstruction &&
        report.finite_schur_inequality_verified &&
        report.finite_half_gap_decay_inequality_verified;
    return report;
}

}  // namespace lemma
