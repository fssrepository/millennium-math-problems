#include "local_quartic_closure_reporter.hpp"

#include "state_analysis.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace lemma {
namespace {

std::string state_filename(const std::string& directory, int cutoff) {
    return (std::filesystem::path(directory) /
            ("K" + std::to_string(cutoff) + ".tsv")).string();
}

std::string objective_formula(const std::string& objective) {
    if (objective == "closure-ratio") {
        return "maximize |K+G| E^(1/4) / (Z^(7/4) P)";
    }
    if (objective == "lqc3-ratio") {
        return "maximize |K+G| / (Z^(5/4) P^(3/4))";
    }
    if (objective == "signed-lqc3-ratio") {
        return "maximize (K+G) / (Z^(5/4) P^(3/4))";
    }
    if (objective == "remainder-envelope-ratio") {
        return "maximize the exact double-square upper envelope / (Z^(5/4) P^(3/4))";
    }
    if (objective == "remainder-absorption-ratio") {
        return "maximize [K+G+(1-theta)||A^(1/2)(B-cAu)||^2] / (Z^(5/4) P^(3/4))";
    }
    if (objective == "shape-power-ratio") {
        return "maximize c_selected^2 x_full^(2p), p=shape_power";
    }
    if (objective == "projective-coherence-ratio") {
        return "maximize ||sum_sigma B_sigma||_2^2 / sum_sigma ||B_sigma||_2^2";
    }
    if (objective == "projective-stretching-ratio") {
        return "maximize |<Au,sum_sigma B_sigma>|^2 / (||Au||_2^2 sum_sigma ||B_sigma||_2^2)";
    }
    if (objective == "projective-cross-power-ratio") {
        return "maximize |(K_cross+G_cross) S_full|^2 / (Z^4 P^4) after subtracting every same-projective-ray quartet";
    }
    if (objective == "projective-open-power-ratio") {
        return "maximize |(K_open+G_open) S_full|^2 / (Z^4 P^4) after subtracting the complete quartet internal to a fixed primitive-shape-height core";
    }
    if (objective == "projective-height-stretching-ratio") {
        return "maximize |<Au,B_(H,2H](u,u)>|^2 / (Z sum_sigma ||A^(1/2)B_sigma(u,u)||_2^2)";
    }
    if (objective == "projective-height-power-ratio") {
        return "maximize |J_(H,2H] S_full|^2 / (Z^4 P^4), where J_(H,2H] is the complete quartet internal to one dyadic primitive-height shell";
    }
    if (objective == "projective-height-outer-power-ratio") {
        return "maximize |S_full|^2 (sum_j ||A^(1/2)B_j(u,u)||_2^2)^2 / (Z^4 P^4) over dyadic primitive-height shells";
    }
    if (objective == "projective-height-envelope-ratio") {
        return "maximize the square of the exact five-component absolute dyadic height-matrix envelope times |S_full| / (Z^2 P^2)";
    }
    if (objective ==
        "projective-height-commutator-envelope-ratio") {
        return "maximize the squared dyadic height envelope after pairing the outer-square and advected terms into their exact commutator";
    }
    if (objective ==
        "projective-height-dynamic-envelope-ratio") {
        return "maximize the squared dyadic height envelope after jointly pairing the outer-square, advected, and nested dynamical terms";
    }
    if (objective ==
        "projective-height-commutator-coercivity-ratio") {
        return "maximize the squared commutator-paired height envelope divided by the squared outer H1 shell weight";
    }
    if (objective ==
        "projective-height-dynamic-coercivity-ratio") {
        return "maximize the squared height envelope after jointly pairing the outer, advected, and nested dynamical terms, divided by the squared outer H1 shell weight";
    }
    if (objective ==
        "projective-palinstrophy-normalization-ratio") {
        return "maximize the squared exact palinstrophy-normalization contribution |S_full 3 S_selected T_selected/(2P)|/(Z^2 P^2)";
    }
    if (objective ==
        "projective-open-palinstrophy-normalization-ratio") {
        return "maximize the squared palinstrophy-normalization contribution remaining after subtracting the complete fixed primitive-height core";
    }
    if (objective ==
        "projective-tail-stretching-alignment-ratio") {
        return "maximize |<Au,B_tail(u,u)>|^2 / (Z ||A^(1/2)B_tail(u,u)||_2^2) above a fixed primitive-height core";
    }
    if (objective ==
        "projective-normalization-alignment-ratio") {
        return "maximize the squared product of the selected-stretching H1 alignment and tail-palinstrophy-cross H2 alignment in the canonical normalization channel";
    }
    if (objective ==
        "projective-normalization-cauchy-ratio") {
        return "maximize the squared exact Cauchy majorant for the canonical selected-stretching times tail-palinstrophy-cross normalization channel";
    }
    if (objective ==
        "projective-selected-stretching-tail-cross-ratio") {
        return "maximize the squared normalized |S_full s_selected t_tail| channel in the two-term open palinstrophy-normalization factorization";
    }
    if (objective ==
        "projective-core-stretching-tail-cross-ratio") {
        return "maximize the squared normalized |S_full s_core t_tail| channel in the open palinstrophy-normalization factorization";
    }
    if (objective ==
        "projective-tail-stretching-core-cross-ratio") {
        return "maximize the squared normalized |S_full s_tail t_core| channel in the open palinstrophy-normalization factorization";
    }
    if (objective ==
        "projective-tail-stretching-tail-cross-ratio") {
        return "maximize the squared normalized |S_full s_tail t_tail| channel in the open palinstrophy-normalization factorization";
    }
    if (objective == "signed-closure-ratio") {
        return "maximize (K+G) E^(1/4) / (Z^(7/4) P)";
    }
    if (objective == "block-ratio") {
        return "maximize selected closed bracket times the full local SLD shape factor";
    }
    if (objective == "mixed-ratio") {
        return "maximize mixed bracket times the full local SLD shape factor";
    }
    if (objective == "terminal-sld-ratio") {
        return "maximize the terminal local SLD ratio with k0 and B0 frozen from the initial state";
    }
    if (objective == "maximum-sld-ratio") {
        return "maximize the largest local SLD ratio on the trajectory with k0 and B0 frozen from the initial state";
    }
    return "maximize 4 S^3 Z P (K+G) / (k0 (S^4 Z^2 P + B0 Z^3 P^4))";
}

void write_json(const LocalQuarticClosureAdversaryReport& report,
                const LocalQuarticClosureAdversaryOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write local quartic closure certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-quartic-closure-adversary-v1\",\n"
        << "  \"search_objective\": \"" << options.objective << "\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"optimized_formula\": \""
        << objective_formula(options.objective) << "\",\n"
        << "  \"closure_objective\": \"maximize |K+G| / (k0 B0^(1/4) Z^(5/4) P^(3/4))\",\n"
        << "  \"direct_sld_objective\": \"maximize 4 S^3 Z P (K+G) / (k0 (S^4 Z^2 P + B0 Z^3 P^4))\",\n"
        << "  \"implemented_static_scale\": \"Z^(7/4) P / E^(1/4)\",\n"
        << "  \"optimizer\": \"" << options.method << "\",\n"
        << "  \"backend\": \"" << options.backend << "\",\n"
        << "  \"initial_profile\": \""
        << options.initial_profile << "\",\n"
        << "  \"warm_state_path\": \"" << options.warm_state_path
        << "\",\n"
        << "  \"lean_diagnostics\": "
        << (options.lean_diagnostics ? "true" : "false") << ",\n"
        << "  \"preserves_sparse_warm_layout\": "
        << (options.preserve_warm_layout ? "true" : "false") << ",\n"
        << "  \"complete_galerkin_cutoff\": "
        << (options.preserve_warm_layout ? "false" : "true") << ",\n"
        << "  \"gradient\": \"exact discrete reverse mode; selected local objectives use direct triads and RK4 follows the requested backend\",\n"
        << "  \"workers\": " << report.workers << ",\n"
        << "  \"restarts_per_cutoff\": " << report.restarts << ",\n"
        << "  \"iterations_per_restart\": " << report.iterations << ",\n"
        << "  \"trajectory_steps\": " << report.trajectory_steps << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(report.viscosity) << ",\n"
        << "  \"time_step\": "
        << static_cast<double>(report.time_step) << ",\n"
        << "  \"absorption_theta\": "
        << static_cast<double>(report.absorption_theta) << ",\n"
        << "  \"shape_power\": " << report.shape_power << ",\n"
        << "  \"projective_core_maximum_height\": "
        << report.projective_core_maximum_height << ",\n"
        << "  \"sobolev_order\": " << report.sobolev_order << ",\n"
        << "  \"sobolev_cap\": "
        << static_cast<double>(report.sobolev_cap) << ",\n"
        << "  \"fitted_cutoff_slope\": "
        << static_cast<double>(report.fitted_cutoff_slope) << ",\n"
        << "  \"maximum_constant_ratio\": "
        << static_cast<double>(report.maximum_constant_ratio) << ",\n"
        << "  \"maximum_objective\": "
        << static_cast<double>(report.maximum_objective) << ",\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_search_is_not_a_proof\": true,\n"
        << "  \"remaining_requirement\": \"derive a cutoff-independent analytic bound for the two-entry bracket\",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        const auto& winner = row.winner;
        const auto& value = winner.value;
        output << "    {\"cutoff\": " << row.cutoff
            << ", \"initial_objective\": "
            << static_cast<double>(winner.initial_objective)
            << ", \"optimized_objective\": "
            << static_cast<double>(winner.objective)
            << ", \"objective_gain\": "
            << static_cast<double>(row.objective_gain)
            << ", \"initial_constant_ratio\": "
            << static_cast<double>(winner.initial_constant_ratio)
            << ", \"optimized_constant_ratio\": "
            << static_cast<double>(value.constant_ratio)
            << ", \"lqc3_target_ratio\": "
            << static_cast<double>(value.lqc3_target_ratio)
            << ", \"signed_lqc3_target_ratio\": "
            << static_cast<double>(value.signed_lqc3_target_ratio)
            << ", \"remainder_upper_envelope_ratio\": "
            << static_cast<double>(winner.remainder_envelope_ratio)
            << ", \"remainder_absorption_ratio\": "
            << static_cast<double>(winner.remainder_absorption_ratio)
            << ", \"shape_power_absolute_product\": "
            << static_cast<double>(winner.shape_power_absolute_product)
            << ", \"shape_power_normalized_stretching\": "
            << static_cast<double>(
                   winner.shape_power_normalized_stretching)
            << ", \"projective_coherence_ratio\": "
            << static_cast<double>(winner.projective_coherence_ratio)
            << ", \"projective_coherence_amplification\": "
            << static_cast<double>(
                   winner.projective_coherence_amplification)
            << ", \"projective_coherence_shape_count\": "
            << winner.projective_coherence_shape_count
            << ", \"projective_stretching_ratio\": "
            << static_cast<double>(winner.projective_stretching_ratio)
            << ", \"projective_stretching_alignment_squared\": "
            << static_cast<double>(
                   winner.projective_stretching_alignment_squared)
            << ", \"projective_stretching_reconstruction_error\": "
            << static_cast<double>(
                   winner.projective_stretching_reconstruction_error)
            << ", \"projective_tail_stretching_alignment_squared\": "
            << static_cast<double>(
                   winner.projective_tail_stretching_alignment_squared)
            << ", \"projective_cross_power_absolute\": "
            << static_cast<double>(winner.projective_cross_power_absolute)
            << ", \"projective_cross_bracket\": "
            << static_cast<double>(winner.projective_cross_bracket)
            << ", \"projective_diagonal_bracket\": "
            << static_cast<double>(winner.projective_diagonal_bracket)
            << ", \"projective_open_power_absolute\": "
            << static_cast<double>(
                   winner.projective_open_power_absolute)
            << ", \"projective_open_bracket\": "
            << static_cast<double>(winner.projective_open_bracket)
            << ", \"projective_fixed_core_bracket\": "
            << static_cast<double>(
                   winner.projective_fixed_core_bracket)
            << ", \"projective_height_stretching_ratio\": "
            << static_cast<double>(
                   winner.projective_height_stretching_ratio)
            << ", \"projective_height_h1_synthesis_ratio\": "
            << static_cast<double>(
                   winner.projective_height_h1_synthesis_ratio)
            << ", \"projective_height_stretching_alignment_squared\": "
            << static_cast<double>(
                   winner.projective_height_stretching_alignment_squared)
            << ", \"projective_height_shape_count\": "
            << winner.projective_height_shape_count
            << ", \"projective_height_power_absolute\": "
            << static_cast<double>(
                   winner.projective_height_power_absolute)
            << ", \"projective_height_internal_bracket\": "
            << static_cast<double>(
                   winner.projective_height_internal_bracket)
            << ", \"projective_height_outer_power_absolute\": "
            << static_cast<double>(
                   winner.projective_height_outer_power_absolute)
            << ", \"projective_height_outer_h1_sum\": "
            << static_cast<double>(
                   winner.projective_height_outer_h1_sum)
            << ", \"projective_height_component_envelope_absolute\": "
            << static_cast<double>(
                   winner.projective_height_component_envelope_absolute)
            << ", \"projective_height_commutator_envelope_absolute\": "
            << static_cast<double>(
                   winner.projective_height_commutator_envelope_absolute)
            << ", \"projective_height_dynamic_envelope_absolute\": "
            << static_cast<double>(
                   winner.projective_height_dynamic_envelope_absolute)
            << ", \"projective_height_commutator_coercivity_ratio\": "
            << static_cast<double>(
                   winner.projective_height_commutator_coercivity_ratio)
            << ", \"projective_height_dynamic_coercivity_ratio\": "
            << static_cast<double>(
                   winner.projective_height_dynamic_coercivity_ratio)
            << ", \"projective_palinstrophy_normalization_power_one\": "
            << static_cast<double>(
                   winner
                       .projective_palinstrophy_normalization_power_one)
            << ", \"projective_height_component_bracket_envelope\": "
            << static_cast<double>(
                   winner.projective_height_component_bracket_envelope)
            << ", \"projective_height_pair_count\": "
            << winner.projective_height_pair_count
            << ", \"projective_height_active_shell_count\": "
            << winner.projective_height_active_shell_count
            << ", \"squared_lqc3_target_ratio\": "
            << static_cast<double>(value.squared_lqc3_target_ratio)
            << ", \"signed_constant_ratio\": "
            << static_cast<double>(value.signed_constant_ratio)
            << ", \"normalized_stretching_ratio\": "
            << static_cast<double>(value.normalized_stretching_ratio)
            << ", \"signed_shape_factor\": "
            << static_cast<double>(value.signed_shape_factor)
            << ", \"factorized_local_sld_ratio\": "
            << static_cast<double>(value.factorized_local_sld_ratio)
            << ", \"factorization_relative_error\": "
            << static_cast<double>(value.factorization_relative_error)
            << ", \"common_block_objective\": "
            << (winner.common_block_objective ? "true" : "false")
            << ", \"common_block_constant_ratio\": "
            << static_cast<double>(
                   winner.common_block_value.block_constant_ratio)
            << ", \"common_block_shape_factor\": "
            << static_cast<double>(
                   winner.common_block_value.common_shape_factor)
            << ", \"common_block_sld_ratio\": "
            << static_cast<double>(
                   winner.common_block_value.block_sld_ratio)
            << ", \"common_block_reconstruction_error\": "
            << static_cast<double>(
                   winner.common_block_value.ratio_reconstruction_error)
            << ", \"refined_objective\": "
            << static_cast<double>(winner.refined_objective)
            << ", \"time_step_relative_error\": "
            << static_cast<double>(winner.time_step_relative_error)
            << ", \"frozen_initial_frequency\": "
            << static_cast<double>(winner.frozen_initial_frequency)
            << ", \"frozen_initial_ep_shift\": "
            << static_cast<double>(winner.frozen_initial_ep_shift)
            << ", \"objective_step\": " << winner.objective_step
            << ", \"refined_objective_step\": "
            << winner.refined_objective_step
            << ", \"improvement_factor\": "
            << static_cast<double>(row.improvement_factor)
            << ", \"signed_two_entry_bracket\": "
            << static_cast<double>(value.signed_two_entry_bracket)
            << ", \"candidate_scale\": "
            << static_cast<double>(value.candidate_scale)
            << ", \"signed_stretching\": "
            << static_cast<double>(value.signed_stretching)
            << ", \"local_polynomial_numerator\": "
            << static_cast<double>(value.local_polynomial_numerator)
            << ", \"local_polynomial_denominator\": "
            << static_cast<double>(value.local_polynomial_denominator)
            << ", \"signed_local_sld_ratio\": "
            << static_cast<double>(value.signed_local_sld_ratio)
            << ", \"energy\": " << static_cast<double>(value.energy)
            << ", \"enstrophy\": "
            << static_cast<double>(value.enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(value.palinstrophy)
            << ", \"warm_lift_constant_ratio\": "
            << static_cast<double>(row.warm_lift_constant_ratio)
            << ", \"warm_lift_objective\": "
            << static_cast<double>(row.warm_lift_objective)
            << ", \"projection_residual\": "
            << static_cast<double>(row.projection_residual)
            << ", \"projected_gradient_norm\": "
            << static_cast<double>(winner.final_projected_gradient_norm)
            << ", \"sobolev_value\": "
            << static_cast<double>(winner.sobolev_value)
            << ", \"winning_restart\": " << winner.restart
            << ", \"winning_seed\": \"" << winner.seed
            << "\", \"warm_continuation_won\": "
            << (winner.warm_continuation ? "true" : "false")
            << ", \"accepted_steps\": " << winner.accepted_steps
            << ", \"evaluations\": " << winner.evaluations
            << ", \"state_path\": \"" << row.state_path
            << "\", \"restart_constant_ratios\": [";
        for (std::size_t restart = 0;
             restart < row.restart_constant_ratios.size(); ++restart) {
            output << static_cast<double>(
                row.restart_constant_ratios[restart]);
            if (restart + 1 != row.restart_constant_ratios.size()) {
                output << ", ";
            }
        }
        output << "], \"restart_objectives\": [";
        for (std::size_t restart = 0;
             restart < row.restart_objectives.size(); ++restart) {
            output << static_cast<double>(row.restart_objectives[restart]);
            if (restart + 1 != row.restart_objectives.size()) {
                output << ", ";
            }
        }
        output << "]}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
}

}  // namespace

void LocalQuarticClosureReporter::write_artifacts(
    LocalQuarticClosureAdversaryReport& report,
    const LocalQuarticClosureAdversaryOptions& options) {
    if (options.certificate_path.empty() ||
        options.state_directory.empty()) {
        throw std::invalid_argument(
            "closure adversary requires certificate and state directory");
    }
    std::filesystem::create_directories(options.state_directory);
    for (auto& row : report.rows) {
        row.state_path = state_filename(
            options.state_directory, row.cutoff);
        std::ostringstream metadata;
        metadata << "exact-gradient local quartic closure adversary; "
                 << "cutoff=" << row.cutoff
                 << "; objective=" << options.objective
                 << "; selection=" << options.selection
                 << "; backend=" << options.backend
                 << "; objective_value="
                 << static_cast<double>(row.winner.objective)
                 << "; constant_ratio="
                 << static_cast<double>(
                        row.winner.value.constant_ratio)
                 << "; candidate_lemma_proved=false";
        SpectralStateWriter::write_tsv(
            row.state_path, row.winner.state, metadata.str());
    }
    write_json(report, options);
}

void LocalQuarticClosureReporter::print_summary(
    const LocalQuarticClosureAdversaryReport& report,
    std::ostream& out) {
    out << "local quartic closure exact-gradient adversary"
        << " objective=" << report.objective
        << " backend=" << report.backend
        << " profile=" << report.initial_profile
        << " workers=" << report.workers
        << " restarts=" << report.restarts
        << " iterations=" << report.iterations << '\n'
        << "absorption_theta="
        << static_cast<double>(report.absorption_theta) << '\n'
        << "cutoff,initial_objective,optimized_objective,gain,"
           "closure_C,lqc3_C,signed_lqc3_C,envelope_C,absorption_C,projective_coherence,projective_amplification,projective_stretching,projective_alignment_squared,projective_cross_power,signed_S,warm_lift_objective,projection_residual,"
           "projective_open_power,gradient_norm,time_refinement_error,accepted,evaluations,seed\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ','
            << static_cast<double>(row.winner.initial_objective) << ','
            << static_cast<double>(row.winner.objective) << ','
            << static_cast<double>(row.objective_gain) << ','
            << static_cast<double>(row.winner.value.constant_ratio) << ','
            << static_cast<double>(row.winner.value.lqc3_target_ratio) << ','
            << static_cast<double>(
                   row.winner.value.signed_lqc3_target_ratio) << ','
            << static_cast<double>(
                   row.winner.remainder_envelope_ratio) << ','
            << static_cast<double>(
                   row.winner.remainder_absorption_ratio) << ','
            << static_cast<double>(
                   row.winner.projective_coherence_ratio) << ','
            << static_cast<double>(
                   row.winner.projective_coherence_amplification) << ','
            << static_cast<double>(
                   row.winner.projective_stretching_ratio) << ','
            << static_cast<double>(
                   row.winner.projective_stretching_alignment_squared) << ','
            << static_cast<double>(
                   row.winner.projective_cross_power_absolute) << ','
            << static_cast<double>(row.winner.value.signed_stretching) << ','
            << static_cast<double>(row.warm_lift_objective) << ','
            << static_cast<double>(row.projection_residual) << ','
            << static_cast<double>(
                   row.winner.projective_open_power_absolute) << ','
            << static_cast<double>(
                   row.winner.final_projected_gradient_norm) << ','
            << static_cast<double>(
                   row.winner.time_step_relative_error) << ','
            << row.winner.accepted_steps << ','
            << row.winner.evaluations << ','
            << row.winner.seed << '\n';
    }
    out << "fitted cutoff slope="
        << static_cast<double>(report.fitted_cutoff_slope)
        << " maximum C="
        << static_cast<double>(report.maximum_constant_ratio)
        << " maximum objective="
        << static_cast<double>(report.maximum_objective)
        << " candidate lemma proved=false\n";
}

}  // namespace lemma
