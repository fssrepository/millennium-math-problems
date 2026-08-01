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
        << "  \"closure_objective\": \"maximize |K+G| / (k0 B0^(1/4) Z^(5/4) P^(3/4))\",\n"
        << "  \"direct_sld_objective\": \"maximize 4 S^3 Z P (K+G) / (k0 (S^4 Z^2 P + B0 Z^3 P^4))\",\n"
        << "  \"implemented_static_scale\": \"Z^(7/4) P / E^(1/4)\",\n"
        << "  \"optimizer\": \"" << options.method << "\",\n"
        << "  \"warm_state_path\": \"" << options.warm_state_path
        << "\",\n"
        << "  \"gradient\": \"exact discrete reverse mode through local Galerkin triads\",\n"
        << "  \"workers\": " << report.workers << ",\n"
        << "  \"restarts_per_cutoff\": " << report.restarts << ",\n"
        << "  \"iterations_per_restart\": " << report.iterations << ",\n"
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
        << " workers=" << report.workers
        << " restarts=" << report.restarts
        << " iterations=" << report.iterations << '\n'
        << "cutoff,initial_objective,optimized_objective,gain,"
           "closure_C,signed_S,warm_lift_objective,projection_residual,"
           "gradient_norm,accepted,evaluations,seed\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ','
            << static_cast<double>(row.winner.initial_objective) << ','
            << static_cast<double>(row.winner.objective) << ','
            << static_cast<double>(row.objective_gain) << ','
            << static_cast<double>(row.winner.value.constant_ratio) << ','
            << static_cast<double>(row.winner.value.signed_stretching) << ','
            << static_cast<double>(row.warm_lift_objective) << ','
            << static_cast<double>(row.projection_residual) << ','
            << static_cast<double>(
                   row.winner.final_projected_gradient_norm) << ','
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
