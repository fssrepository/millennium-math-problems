#include "local_sld_projective_height_matrix.hpp"

#include "local_sld_projective_height_tail_summary.hpp"
#include "local_sld_projective_height_schur_summary.hpp"
#include "projective_advection_decomposition.hpp"
#include "projective_height_shell_partition.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

struct HeightShellData {
    LocalSldProjectiveHeightShellSummary summary;
    std::vector<std::size_t> group_indices;
    std::vector<std::size_t> target_multiplicity;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralIncrement transported_au;
    SpectralIncrement commutator;
};

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective height-matrix Laplacian layout mismatch");
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

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective height-matrix pairing layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal hminus1_norm_squared(
    const SpectralState& state,
    const SpectralIncrement& field) {
    if (field.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective height-matrix H-1 layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < field.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        if (wave2 > 0.0L) {
            result += std::real(dot_hermitian(
                field[mode], field[mode])) / wave2;
        }
    }
    return result;
}

SpectralReal relative_error(
    SpectralReal computed,
    SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

TriadSelection selection_for(
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (exclude_triple_family) {
        return exclude_signature_123
            ? TriadSelection::
                  local_without_equal_low_double_triple_and_signature(
                      1, 2, 3)
            : TriadSelection::local_without_equal_low_double_triple();
    }
    return exclude_signature_123
        ? TriadSelection::
              local_without_equal_low_doubling_and_signature(1, 2, 3)
        : TriadSelection::local_without_equal_low_doubling();
}

void finalize_entry(
    LocalSldProjectiveHeightMatrixEntry& entry,
    SpectralReal power_one_scale) {
    entry.bracket = entry.outer_square +
        entry.advected_commutator + entry.advecting_nested +
        entry.enstrophy_normalization +
        entry.palinstrophy_normalization;
    entry.power_one = entry.bracket * power_one_scale;
    entry.absolute_component_power_one_envelope =
        (std::abs(entry.outer_square) +
         std::abs(entry.advected_commutator) +
         std::abs(entry.advecting_nested) +
         std::abs(entry.enstrophy_normalization) +
         std::abs(entry.palinstrophy_normalization)) *
        std::abs(power_one_scale);
    entry.commutator_paired_power_one_envelope =
        (std::abs(entry.outer_square + entry.advected_commutator) +
         std::abs(entry.advecting_nested) +
         std::abs(entry.enstrophy_normalization) +
         std::abs(entry.palinstrophy_normalization)) *
        std::abs(power_one_scale);
    entry.dynamic_paired_power_one_envelope =
        (std::abs(entry.outer_square + entry.advected_commutator +
                  entry.advecting_nested) +
         std::abs(entry.enstrophy_normalization) +
         std::abs(entry.palinstrophy_normalization)) *
        std::abs(power_one_scale);
}

void write_json(
    const LocalSldProjectiveHeightMatrixReport& report,
    const LocalSldProjectiveHeightTailReport& tail,
    const LocalSldProjectiveHeightSchurReport& schur,
    const LocalSldProjectiveHeightMatrixCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-matrix-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"excludes_signature_123\": "
        << (report.excludes_signature_123 ? "true" : "false")
        << ",\n"
        << "  \"excludes_triple_family\": "
        << (report.excludes_triple_family ? "true" : "false")
        << ",\n"
        << "  \"power_one_scale\": "
        << static_cast<double>(report.power_one_scale) << ",\n"
        << "  \"selected_power_one\": "
        << static_cast<double>(report.selected_power_one) << ",\n"
        << "  \"reconstructed_power_one\": "
        << static_cast<double>(report.reconstructed_power_one) << ",\n"
        << "  \"absolute_power_one_sum\": "
        << static_cast<double>(report.absolute_power_one_sum) << ",\n"
        << "  \"effective_height_pairs\": "
        << static_cast<double>(report.effective_height_pairs) << ",\n"
        << "  \"dominant_height_pair_fraction\": "
        << static_cast<double>(report.dominant_height_pair_fraction)
        << ",\n"
        << "  \"signed_height_pair_alignment\": "
        << static_cast<double>(report.signed_height_pair_alignment)
        << ",\n"
        << "  \"shells\": [\n";
    for (std::size_t index = 0; index < report.shells.size(); ++index) {
        const auto& shell = report.shells[index];
        output << "    {\"shell\": " << shell.shell
            << ", \"minimum_height\": " << shell.minimum_height
            << ", \"maximum_height\": " << shell.maximum_height
            << ", \"shape_count\": " << shell.shape_count
            << ", \"interaction_count\": "
            << shell.interaction_count
            << ", \"target_mode_count\": "
            << shell.target_mode_count
            << ", \"target_multiplicity_l2_squared\": "
            << static_cast<double>(
                   shell.target_multiplicity_l2_squared)
            << ", \"stretching\": "
            << static_cast<double>(shell.stretching)
            << ", \"palinstrophy_cross\": "
            << static_cast<double>(shell.palinstrophy_cross)
            << ", \"aggregate_l2_norm2\": "
            << static_cast<double>(shell.aggregate_l2_norm2)
            << ", \"aggregate_h1_norm2\": "
            << static_cast<double>(shell.aggregate_h1_norm2)
            << ", \"aggregate_h2_norm2\": "
            << static_cast<double>(shell.aggregate_h2_norm2)
            << ", \"square_function_l2_norm2\": "
            << static_cast<double>(shell.square_function_l2_norm2)
            << ", \"square_function_h1_norm2\": "
            << static_cast<double>(shell.square_function_h1_norm2)
            << ", \"square_function_h2_norm2\": "
            << static_cast<double>(shell.square_function_h2_norm2)
            << ", \"l2_synthesis_ratio\": "
            << static_cast<double>(shell.l2_synthesis_ratio)
            << ", \"h1_synthesis_ratio\": "
            << static_cast<double>(shell.h1_synthesis_ratio)
            << ", \"h2_synthesis_ratio\": "
            << static_cast<double>(shell.h2_synthesis_ratio)
            << ", \"stretching_alignment_squared\": "
            << static_cast<double>(shell.stretching_alignment_squared)
            << ", \"stretching_h1_alignment_squared\": "
            << static_cast<double>(
                   shell.stretching_h1_alignment_squared)
            << ", \"h1_synthesis_stretching_product\": "
            << static_cast<double>(
                   shell.h1_synthesis_stretching_product)
            << ", \"palinstrophy_cross_alignment_squared\": "
            << static_cast<double>(
                   shell.palinstrophy_cross_alignment_squared)
            << ", \"commutator_hminus1_norm2\": "
            << static_cast<double>(
                   shell.commutator_hminus1_norm2)
            << ", \"commutator_to_outer_ratio_squared\": "
            << static_cast<double>(
                   shell.commutator_to_outer_ratio_squared)
            << ", \"diagonal_commutator_alignment_squared\": "
            << static_cast<double>(
                   shell.diagonal_commutator_alignment_squared)
            << '}'
            << (index + 1 == report.shells.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"height_pair_entries\": [\n";
    for (std::size_t index = 0; index < report.entries.size(); ++index) {
        const auto& entry = report.entries[index];
        output << "    {\"first_shell\": " << entry.first_shell
            << ", \"second_shell\": " << entry.second_shell
            << ", \"first_height_range\": ["
            << entry.first_minimum_height << ", "
            << entry.first_maximum_height << ']'
            << ", \"second_height_range\": ["
            << entry.second_minimum_height << ", "
            << entry.second_maximum_height << ']'
            << ", \"outer_square\": "
            << static_cast<double>(entry.outer_square)
            << ", \"advected_commutator\": "
            << static_cast<double>(entry.advected_commutator)
            << ", \"advecting_nested\": "
            << static_cast<double>(entry.advecting_nested)
            << ", \"enstrophy_normalization\": "
            << static_cast<double>(entry.enstrophy_normalization)
            << ", \"palinstrophy_normalization\": "
            << static_cast<double>(entry.palinstrophy_normalization)
            << ", \"power_one\": "
            << static_cast<double>(entry.power_one)
            << ", \"absolute_component_power_one_envelope\": "
            << static_cast<double>(
                   entry.absolute_component_power_one_envelope)
            << ", \"commutator_paired_power_one_envelope\": "
            << static_cast<double>(
                   entry.commutator_paired_power_one_envelope)
            << ", \"dynamic_paired_power_one_envelope\": "
            << static_cast<double>(
                   entry.dynamic_paired_power_one_envelope)
            << ", \"shared_target_mode_count\": "
            << entry.shared_target_mode_count
            << ", \"target_incidence_cosine\": "
            << static_cast<double>(entry.target_incidence_cosine)
            << ", \"absolute_fraction\": "
            << static_cast<double>(entry.absolute_fraction) << '}'
            << (index + 1 == report.entries.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"cumulative_height_cuts\": [\n";
    for (std::size_t index = 0; index < tail.rows.size(); ++index) {
        const auto& row = tail.rows[index];
        output << "    {\"core_maximum_height\": "
            << row.core_maximum_height
            << ", \"core_internal_power_one\": "
            << static_cast<double>(row.core_internal_power_one)
            << ", \"core_tail_power_one\": "
            << static_cast<double>(row.core_tail_power_one)
            << ", \"tail_internal_power_one\": "
            << static_cast<double>(row.tail_internal_power_one)
            << ", \"open_power_one\": "
            << static_cast<double>(row.open_power_one)
            << ", \"open_outer_square_power_one\": "
            << static_cast<double>(row.open_outer_square_power_one)
            << ", \"open_advected_commutator_power_one\": "
            << static_cast<double>(
                   row.open_advected_commutator_power_one)
            << ", \"open_advecting_nested_power_one\": "
            << static_cast<double>(row.open_advecting_nested_power_one)
            << ", \"open_enstrophy_normalization_power_one\": "
            << static_cast<double>(
                   row.open_enstrophy_normalization_power_one)
            << ", \"open_palinstrophy_normalization_power_one\": "
            << static_cast<double>(
                   row.open_palinstrophy_normalization_power_one)
            << ", \"open_absolute_power_one_sum\": "
            << static_cast<double>(row.open_absolute_power_one_sum)
            << ", \"open_effective_height_pairs\": "
            << static_cast<double>(row.open_effective_height_pairs)
            << ", \"dominant_open_pair_fraction\": "
            << static_cast<double>(row.dominant_open_pair_fraction)
            << ", \"open_signed_alignment\": "
            << static_cast<double>(row.open_signed_alignment)
            << ", \"reconstruction_error\": "
            << static_cast<double>(row.reconstruction_error)
            << ", \"component_reconstruction_error\": "
            << static_cast<double>(row.component_reconstruction_error)
            << '}'
            << (index + 1 == tail.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"height_schur_gaps\": [\n";
    for (std::size_t index = 0; index < schur.gaps.size(); ++index) {
        const auto& gap = schur.gaps[index];
        output << "    {\"shell_gap\": " << gap.shell_gap
            << ", \"pair_count\": " << gap.pair_count
            << ", \"unscaled_pair_count\": "
            << gap.unscaled_pair_count
            << ", \"signed_power_one_sum\": "
            << static_cast<double>(gap.signed_power_one_sum)
            << ", \"absolute_power_one_sum\": "
            << static_cast<double>(gap.absolute_power_one_sum)
            << ", \"absolute_component_envelope_sum\": "
            << static_cast<double>(
                   gap.absolute_component_envelope_sum)
            << ", \"maximum_symmetric_geometric_ratio\": "
            << static_cast<double>(
                   gap.maximum_symmetric_geometric_ratio)
            << ", \"commutator_paired_envelope_sum\": "
            << static_cast<double>(
                   gap.commutator_paired_envelope_sum)
            << ", \"commutator_paired_outer_maximum_symmetric_geometric_ratio\": "
            << static_cast<double>(
                   gap
                       .commutator_paired_outer_maximum_symmetric_geometric_ratio)
            << ", \"commutator_term_envelope_sum\": "
            << static_cast<double>(gap.commutator_term_envelope_sum)
            << ", \"remainder_terms_envelope_sum\": "
            << static_cast<double>(gap.remainder_terms_envelope_sum)
            << ", \"commutator_outer_maximum_symmetric_geometric_ratio\": "
            << static_cast<double>(
                   gap
                       .commutator_outer_maximum_symmetric_geometric_ratio)
            << ", \"remainder_outer_maximum_symmetric_geometric_ratio\": "
            << static_cast<double>(
                   gap.remainder_outer_maximum_symmetric_geometric_ratio)
            << ", \"commutator_paired_outer_unscaled_pair_count\": "
            << gap.commutator_paired_outer_unscaled_pair_count << '}'
            << (index + 1 == schur.gaps.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"component_schur_diagnostics\": [\n";
    for (std::size_t index = 0; index < schur.components.size(); ++index) {
        const auto& component = schur.components[index];
        output << "    {\"component\": \"" << component.component
            << "\", \"total_power_one_envelope\": "
            << static_cast<double>(
                   component.total_power_one_envelope)
            << ", \"diagonal_power_one_envelope\": "
            << static_cast<double>(
                   component.diagonal_power_one_envelope)
            << ", \"maximum_off_diagonal_geometric_ratio\": "
            << static_cast<double>(
                   component.maximum_off_diagonal_geometric_ratio)
            << ", \"maximum_weighted_row_sum\": "
            << static_cast<double>(component.maximum_weighted_row_sum)
            << ", \"unscaled_off_diagonal_pair_count\": "
            << component.unscaled_off_diagonal_pair_count
            << ", \"maximum_outer_weighted_diagonal_ratio\": "
            << static_cast<double>(
                   component.maximum_outer_weighted_diagonal_ratio)
            << ", \"maximum_outer_weighted_off_diagonal_ratio\": "
            << static_cast<double>(
                   component.maximum_outer_weighted_off_diagonal_ratio)
            << ", \"maximum_outer_weighted_row_sum\": "
            << static_cast<double>(
                   component.maximum_outer_weighted_row_sum)
            << ", \"outer_unscaled_pair_count\": "
            << component.outer_unscaled_pair_count << '}'
            << (index + 1 == schur.components.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"maximum_symmetric_geometric_ratio\": "
        << static_cast<double>(
               schur.maximum_symmetric_geometric_ratio)
        << ",\n"
        << "  \"maximum_weighted_schur_row_sum\": "
        << static_cast<double>(schur.maximum_weighted_row_sum)
        << ",\n"
        << "  \"total_component_power_one_envelope\": "
        << static_cast<double>(schur.total_component_envelope)
        << ",\n"
        << "  \"diagonal_component_power_one_envelope\": "
        << static_cast<double>(schur.diagonal_component_envelope)
        << ",\n"
        << "  \"weighted_schur_upper_bound\": "
        << static_cast<double>(schur.weighted_schur_upper_bound)
        << ",\n"
        << "  \"weighted_schur_upper_bound_ratio\": "
        << static_cast<double>(schur.upper_bound_ratio)
        << ",\n"
        << "  \"commutator_paired_total_power_one_envelope\": "
        << static_cast<double>(
               schur.commutator_paired_total_envelope)
        << ",\n"
        << "  \"commutator_paired_diagonal_power_one_envelope\": "
        << static_cast<double>(
               schur.commutator_paired_diagonal_envelope)
        << ",\n"
        << "  \"commutator_paired_maximum_weighted_schur_row_sum\": "
        << static_cast<double>(
               schur.commutator_paired_maximum_weighted_row_sum)
        << ",\n"
        << "  \"commutator_paired_weighted_schur_upper_bound\": "
        << static_cast<double>(
               schur.commutator_paired_weighted_schur_upper_bound)
        << ",\n"
        << "  \"commutator_paired_weighted_schur_upper_bound_ratio\": "
        << static_cast<double>(
               schur.commutator_paired_upper_bound_ratio)
        << ",\n"
        << "  \"commutator_paired_unscaled_off_diagonal_pair_count\": "
        << schur.commutator_paired_unscaled_off_diagonal_pair_count
        << ",\n"
        << "  \"finite_commutator_paired_schur_inequality_verified\": "
        << (schur.finite_commutator_paired_schur_inequality_verified
                ? "true" : "false")
        << ",\n"
        << "  \"commutator_paired_outer_power_one_weight\": "
        << static_cast<double>(
               schur.commutator_paired_outer_weight)
        << ",\n"
        << "  \"commutator_paired_outer_maximum_weighted_schur_row_sum\": "
        << static_cast<double>(
               schur
                   .commutator_paired_outer_maximum_weighted_row_sum)
        << ",\n"
        << "  \"commutator_paired_outer_weighted_schur_upper_bound\": "
        << static_cast<double>(
               schur
                   .commutator_paired_outer_weighted_schur_upper_bound)
        << ",\n"
        << "  \"commutator_paired_outer_weighted_schur_upper_bound_ratio\": "
        << static_cast<double>(
               schur.commutator_paired_outer_upper_bound_ratio)
        << ",\n"
        << "  \"commutator_paired_outer_unscaled_off_diagonal_pair_count\": "
        << schur
               .commutator_paired_outer_unscaled_off_diagonal_pair_count
        << ",\n"
        << "  \"finite_commutator_paired_outer_schur_inequality_verified\": "
        << (schur
                    .finite_commutator_paired_outer_schur_inequality_verified
                ? "true" : "false")
        << ",\n"
        << "  \"dynamic_paired_total_power_one_envelope\": "
        << static_cast<double>(schur.dynamic_paired_total_envelope)
        << ",\n"
        << "  \"dynamic_paired_outer_maximum_weighted_schur_row_sum\": "
        << static_cast<double>(
               schur.dynamic_paired_outer_maximum_weighted_row_sum)
        << ",\n"
        << "  \"dynamic_paired_outer_weighted_schur_upper_bound\": "
        << static_cast<double>(
               schur.dynamic_paired_outer_weighted_schur_upper_bound)
        << ",\n"
        << "  \"dynamic_paired_outer_weighted_schur_upper_bound_ratio\": "
        << static_cast<double>(
               schur.dynamic_paired_outer_upper_bound_ratio)
        << ",\n"
        << "  \"dynamic_paired_outer_unscaled_off_diagonal_pair_count\": "
        << schur.dynamic_paired_outer_unscaled_off_diagonal_pair_count
        << ",\n"
        << "  \"finite_dynamic_paired_outer_schur_inequality_verified\": "
        << (schur.finite_dynamic_paired_outer_schur_inequality_verified
                ? "true" : "false")
        << ",\n"
        << "  \"finite_weighted_schur_inequality_verified\": "
        << (schur.finite_schur_inequality_verified
                ? "true" : "false")
        << ",\n"
        << "  \"unscaled_off_diagonal_pair_count\": "
        << schur.unscaled_off_diagonal_pair_count << ",\n"
        << "  \"cutoff_uniform_weighted_schur_bound_proved\": false,\n"
        << "  \"maximum_cumulative_reconstruction_error\": "
        << static_cast<double>(tail.maximum_reconstruction_error)
        << ",\n"
        << "  \"maximum_cumulative_component_reconstruction_error\": "
        << static_cast<double>(
               tail.maximum_component_reconstruction_error)
        << ",\n"
        << "  \"exact_cumulative_height_decomposition\": "
        << (tail.exact_cumulative_decomposition ? "true" : "false")
        << ",\n"
        << "  \"bracket_reconstruction_error\": "
        << static_cast<double>(report.bracket_reconstruction_error)
        << ",\n"
        << "  \"exact_height_matrix_decomposition\": "
        << (report.exact_height_matrix_decomposition
                ? "true" : "false") << ",\n"
        << "  \"uniform_height_matrix_bound_proved\": false,\n"
        << "  \"finite_height_matrix_is_not_a_proof\": true,\n"
        << "  \"remaining_requirement\": \"derive a cutoff-uniform weighted Schur or signed summability estimate for this dyadic primitive-height matrix\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveHeightMatrixReport
LocalSldProjectiveHeightMatrix::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height-matrix threads must be 1..256");
    }
    const TriadSelection selection = selection_for(
        exclude_signature_123, exclude_triple_family);
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    const auto partition = ProjectiveHeightShellPartition::build(groups);
    std::vector<HeightShellData> shells(partition.size());
    for (std::size_t index = 0; index < partition.size(); ++index) {
        HeightShellData& data = shells[index];
        data.summary.shell = partition[index].shell;
        data.summary.minimum_height = partition[index].minimum_height;
        data.summary.maximum_height = partition[index].maximum_height;
        data.summary.shape_count =
            partition[index].group_indices.size();
        data.summary.interaction_count =
            partition[index].interaction_count;
        data.group_indices = partition[index].group_indices;
        data.target_multiplicity.assign(state.waves.size(), 0);
        for (const std::size_t group_index : data.group_indices) {
            for (const InteractionIndex interaction :
                 groups[group_index].interactions) {
                ++data.target_multiplicity[interaction[2]];
            }
        }
        for (const std::size_t multiplicity :
             data.target_multiplicity) {
            if (multiplicity > 0) {
                ++data.summary.target_mode_count;
                data.summary.target_multiplicity_l2_squared +=
                    static_cast<SpectralReal>(multiplicity) *
                    static_cast<SpectralReal>(multiplicity);
            }
        }
    }
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(dynamics, selection).evaluate(state);
    const LocalQuarticClosureObjectiveValue full_local =
        LocalQuarticClosureObjective(
            dynamics, TriadPartition::local).evaluate(state);
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L)) {
        throw std::invalid_argument(
            "projective height matrix requires positive Z and P");
    }
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    auto evaluate = [&](const HeightShellData& shell,
                        const SpectralIncrement& advecting,
                        const SpectralIncrement& advected) {
        return ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, shell.group_indices,
            advecting, advected, threads);
    };
    for (HeightShellData& shell : shells) {
        shell.b = evaluate(shell, state.velocity, state.velocity);
        shell.ab = laplacian_weight(state, shell.b);
        shell.transported_au = evaluate(shell, state.velocity, au);
        shell.commutator = shell.transported_au;
        for (std::size_t mode = 0;
             mode < shell.commutator.size(); ++mode) {
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                shell.commutator[mode][coordinate] -=
                    shell.ab[mode][coordinate];
            }
        }
        shell.summary.stretching = pairing(au, shell.b);
        shell.summary.palinstrophy_cross = pairing(shell.ab, au);
        shell.summary.aggregate_l2_norm2 = pairing(shell.b, shell.b);
        shell.summary.aggregate_h1_norm2 = pairing(shell.b, shell.ab);
        shell.summary.aggregate_h2_norm2 = pairing(shell.ab, shell.ab);
        shell.summary.commutator_hminus1_norm2 =
            hminus1_norm_squared(state, shell.commutator);
        if (shell.summary.aggregate_h1_norm2 > 0.0L) {
            shell.summary.commutator_to_outer_ratio_squared =
                shell.summary.commutator_hminus1_norm2 /
                shell.summary.aggregate_h1_norm2;
        }
        const SpectralReal diagonal_commutator = pairing(
            shell.b, shell.commutator);
        const SpectralReal commutator_alignment_denominator =
            shell.summary.aggregate_h1_norm2 *
            shell.summary.commutator_hminus1_norm2;
        if (commutator_alignment_denominator > 0.0L) {
            shell.summary.diagonal_commutator_alignment_squared =
                diagonal_commutator * diagonal_commutator /
                commutator_alignment_denominator;
        }
        const ProjectiveSquareFunctionNorms square_function =
            ProjectiveAdvectionDecomposition::square_function_norms(
                state, groups, shell.group_indices, threads);
        shell.summary.square_function_l2_norm2 =
            square_function.l2_norm2;
        shell.summary.square_function_h1_norm2 =
            square_function.h1_norm2;
        shell.summary.square_function_h2_norm2 =
            square_function.h2_norm2;
        if (square_function.l2_norm2 > 0.0L) {
            shell.summary.l2_synthesis_ratio =
                shell.summary.aggregate_l2_norm2 /
                square_function.l2_norm2;
        }
        if (square_function.h1_norm2 > 0.0L) {
            shell.summary.h1_synthesis_ratio =
                shell.summary.aggregate_h1_norm2 /
                square_function.h1_norm2;
        }
        if (square_function.h2_norm2 > 0.0L) {
            shell.summary.h2_synthesis_ratio =
                shell.summary.aggregate_h2_norm2 /
                square_function.h2_norm2;
        }
        const SpectralReal stretching_denominator =
            full_local.palinstrophy *
            shell.summary.aggregate_l2_norm2;
        if (stretching_denominator > 0.0L) {
            shell.summary.stretching_alignment_squared =
                shell.summary.stretching * shell.summary.stretching /
                stretching_denominator;
        }
        const SpectralReal stretching_h1_denominator =
            full_local.enstrophy *
            shell.summary.aggregate_h1_norm2;
        if (stretching_h1_denominator > 0.0L) {
            shell.summary.stretching_h1_alignment_squared =
                shell.summary.stretching * shell.summary.stretching /
                stretching_h1_denominator;
        }
        shell.summary.h1_synthesis_stretching_product =
            shell.summary.h1_synthesis_ratio *
            shell.summary.stretching_h1_alignment_squared;
        const SpectralReal cross_denominator =
            full_local.palinstrophy *
            shell.summary.aggregate_h2_norm2;
        if (cross_denominator > 0.0L) {
            shell.summary.palinstrophy_cross_alignment_squared =
                shell.summary.palinstrophy_cross *
                shell.summary.palinstrophy_cross /
                cross_denominator;
        }
    }

    LocalSldProjectiveHeightMatrixReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.threads = threads;
    report.excludes_signature_123 = exclude_signature_123;
    report.excludes_triple_family = exclude_triple_family;
    report.selected_bracket = selected.signed_two_entry_bracket;
    report.power_one_scale =
        full_local.signed_stretching /
        (selected.enstrophy * selected.enstrophy *
         selected.palinstrophy * selected.palinstrophy);
    report.selected_power_one =
        report.selected_bracket * report.power_one_scale;
    SpectralReal square_sum = 0.0L;
    for (std::size_t first = 0; first < shells.size(); ++first) {
        for (std::size_t second = first;
             second < shells.size(); ++second) {
            const HeightShellData& left = shells[first];
            const HeightShellData& right = shells[second];
            LocalSldProjectiveHeightMatrixEntry entry;
            entry.first_shell = static_cast<int>(first);
            entry.second_shell = static_cast<int>(second);
            entry.first_minimum_height = left.summary.minimum_height;
            entry.first_maximum_height = left.summary.maximum_height;
            entry.second_minimum_height = right.summary.minimum_height;
            entry.second_maximum_height = right.summary.maximum_height;
            SpectralReal target_dot = 0.0L;
            for (std::size_t mode = 0;
                 mode < state.waves.size(); ++mode) {
                const std::size_t left_multiplicity =
                    left.target_multiplicity[mode];
                const std::size_t right_multiplicity =
                    right.target_multiplicity[mode];
                if (left_multiplicity > 0 &&
                    right_multiplicity > 0) {
                    ++entry.shared_target_mode_count;
                }
                target_dot +=
                    static_cast<SpectralReal>(left_multiplicity) *
                    static_cast<SpectralReal>(right_multiplicity);
            }
            const SpectralReal incidence_denominator = std::sqrt(
                left.summary.target_multiplicity_l2_squared *
                right.summary.target_multiplicity_l2_squared);
            if (incidence_denominator > 0.0L) {
                entry.target_incidence_cosine =
                    target_dot / incidence_denominator;
            }
            if (first == second) {
                const SpectralIncrement nested = evaluate(
                    left, left.b, state.velocity);
                entry.outer_square = -pairing(left.b, left.ab);
                entry.advected_commutator = pairing(
                    left.b, left.transported_au);
                entry.advecting_nested = -pairing(au, nested);
                entry.enstrophy_normalization =
                    left.summary.stretching *
                    left.summary.stretching /
                    (2.0L * selected.enstrophy);
                entry.palinstrophy_normalization =
                    3.0L * left.summary.stretching *
                    left.summary.palinstrophy_cross /
                    (2.0L * selected.palinstrophy);
            } else {
                const SpectralIncrement left_advects_right = evaluate(
                    left, right.b, state.velocity);
                const SpectralIncrement right_advects_left = evaluate(
                    right, left.b, state.velocity);
                entry.outer_square =
                    -pairing(left.b, right.ab) -
                    pairing(right.b, left.ab);
                entry.advected_commutator =
                    pairing(left.b, right.transported_au) +
                    pairing(right.b, left.transported_au);
                entry.advecting_nested =
                    -pairing(au, left_advects_right) -
                    pairing(au, right_advects_left);
                entry.enstrophy_normalization =
                    left.summary.stretching *
                    right.summary.stretching /
                    selected.enstrophy;
                entry.palinstrophy_normalization =
                    3.0L *
                    (left.summary.stretching *
                         right.summary.palinstrophy_cross +
                     right.summary.stretching *
                         left.summary.palinstrophy_cross) /
                    (2.0L * selected.palinstrophy);
            }
            finalize_entry(entry, report.power_one_scale);
            report.reconstructed_bracket += entry.bracket;
            report.reconstructed_power_one += entry.power_one;
            report.absolute_power_one_sum += std::abs(entry.power_one);
            square_sum += entry.power_one * entry.power_one;
            report.entries.push_back(entry);
        }
    }
    for (const HeightShellData& shell : shells) {
        report.shells.push_back(shell.summary);
    }
    if (square_sum > 0.0L) {
        report.effective_height_pairs =
            report.absolute_power_one_sum *
            report.absolute_power_one_sum / square_sum;
    }
    if (report.absolute_power_one_sum > 0.0L) {
        report.signed_height_pair_alignment =
            std::abs(report.reconstructed_power_one) /
            report.absolute_power_one_sum;
        for (auto& entry : report.entries) {
            entry.absolute_fraction = std::abs(entry.power_one) /
                report.absolute_power_one_sum;
            report.dominant_height_pair_fraction = std::max(
                report.dominant_height_pair_fraction,
                entry.absolute_fraction);
        }
    }
    report.bracket_reconstruction_error = relative_error(
        report.reconstructed_bracket, report.selected_bracket);
    report.exact_height_matrix_decomposition =
        report.bracket_reconstruction_error < 1e-13L;
    return report;
}

LocalSldProjectiveHeightMatrixCliOptions
LocalSldProjectiveHeightMatrixCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightMatrixCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--exclude-123") {
            options.exclude_signature_123 = true;
        } else if (name == "--exclude-triple-family") {
            options.exclude_triple_family = true;
        } else {
            throw std::invalid_argument(
                "unknown projective-height-matrix option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "projective-height-matrix requires state, certificate, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveHeightMatrixCli::print_help(
    std::ostream& out) {
    out << "Dyadic primitive-projective-height quartet matrix options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON matrix\n"
        << "  --threads N           parallel projective interaction workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveHeightMatrixCli::run(
    const LocalSldProjectiveHeightMatrixCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldProjectiveHeightMatrixReport report =
        LocalSldProjectiveHeightMatrix::analyze(
            dynamics, state, options.threads,
            options.exclude_signature_123,
            options.exclude_triple_family);
    const LocalSldProjectiveHeightTailReport tail =
        LocalSldProjectiveHeightTailSummary::summarize(report);
    const LocalSldProjectiveHeightSchurReport schur =
        LocalSldProjectiveHeightSchurSummary::summarize(report);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective height-matrix certificate");
    }
    write_json(report, tail, schur, options, certificate);
    out << std::setprecision(12)
        << "projective height matrix cutoff=" << report.cutoff
        << " shells=" << report.shells.size()
        << " effective_pairs="
        << static_cast<double>(report.effective_height_pairs)
        << " dominant="
        << static_cast<double>(report.dominant_height_pair_fraction)
        << " alignment="
        << static_cast<double>(report.signed_height_pair_alignment)
        << " error="
        << static_cast<double>(report.bracket_reconstruction_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_height_matrix_decomposition &&
            tail.exact_cumulative_decomposition &&
            schur.finite_schur_inequality_verified &&
            schur.finite_dynamic_paired_outer_schur_inequality_verified
        ? 0 : 2;
}

}  // namespace lemma
