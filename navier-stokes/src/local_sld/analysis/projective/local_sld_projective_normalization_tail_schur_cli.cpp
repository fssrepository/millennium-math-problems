#include "local_sld_projective_normalization_tail_schur_cli.hpp"

#include "local_sld_projective_normalization_tail_schur_ledger.hpp"
#include "local_sld_triad_selection.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace lemma {
namespace {

using Report = LocalSldProjectiveNormalizationTailSchurReport;
using ReportMember = SpectralReal Report::*;

SpectralReal fitted_height_slope(
    const std::vector<Report>& reports,
    ReportMember member) {
    SpectralReal count = 0.0L;
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    for (const Report& report : reports) {
        const SpectralReal value = report.*member;
        if (report.core_maximum_height < 1 || !(value > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(report.core_maximum_height));
        const SpectralReal y = std::log(value);
        count += 1.0L;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }
    const SpectralReal denominator =
        count * sum_xx - sum_x * sum_x;
    return count >= 2.0L && std::abs(denominator) > 1e-30L
        ? (count * sum_xy - sum_x * sum_y) / denominator
        : 0.0L;
}

SpectralReal maximum_value(
    const std::vector<Report>& reports,
    ReportMember member) {
    SpectralReal result = 0.0L;
    for (const Report& report : reports) {
        result = std::max(result, report.*member);
    }
    return result;
}

void write_report(
    const LocalSldProjectiveNormalizationTailSchurReport& report,
    const std::string& state_path,
    std::ostream& output) {
    output << std::setprecision(18)
        << "    {\"state_path\": \"" << state_path
        << "\", \"cutoff\": " << report.cutoff
        << ", \"core_maximum_height\": "
        << report.core_maximum_height
        << ", \"selected_shape_count\": "
        << report.selected_shape_count
        << ", \"tail_shape_count\": " << report.tail_shape_count
        << ", \"aggregate_tail_h2_norm2\": "
        << static_cast<double>(report.aggregate_tail_h2_norm2)
        << ", \"reconstructed_tail_h2_norm2\": "
        << static_cast<double>(report.reconstructed_tail_h2_norm2)
        << ", \"diagonal_tail_h2_norm2\": "
        << static_cast<double>(report.diagonal_tail_h2_norm2)
        << ", \"absolute_gram_tail_h2\": "
        << static_cast<double>(report.absolute_gram_tail_h2)
        << ", \"maximum_normalized_absolute_gram_row_sum\": "
        << static_cast<double>(
               report.maximum_normalized_absolute_gram_row_sum)
        << ", \"schur_tail_h2_upper_bound\": "
        << static_cast<double>(report.schur_tail_h2_upper_bound)
        << ", \"maximum_half_decay_weighted_correlation\": "
        << static_cast<double>(
               report.maximum_half_decay_weighted_correlation)
        << ", \"half_decay_implied_row_bound\": "
        << static_cast<double>(report.half_decay_implied_row_bound)
        << ", \"half_decay_tail_h2_upper_bound\": "
        << static_cast<double>(report.half_decay_tail_h2_upper_bound)
        << ", \"aggregate_to_diagonal_ratio\": "
        << static_cast<double>(report.aggregate_to_diagonal_ratio)
        << ", \"aggregate_to_absolute_gram_ratio\": "
        << static_cast<double>(report.aggregate_to_absolute_gram_ratio)
        << ", \"aggregate_to_schur_bound_ratio\": "
        << static_cast<double>(report.aggregate_to_schur_bound_ratio)
        << ", \"tail_h2_reconstruction_error\": "
        << static_cast<double>(report.tail_h2_reconstruction_error)
        << ", \"objective_tail_h2_reconstruction_error\": "
        << static_cast<double>(
               report.objective_tail_h2_reconstruction_error)
        << ", \"full_stretching\": "
        << static_cast<double>(report.full_stretching)
        << ", \"enstrophy\": "
        << static_cast<double>(report.enstrophy)
        << ", \"palinstrophy\": "
        << static_cast<double>(report.palinstrophy)
        << ", \"selected_aggregate_h1_norm2\": "
        << static_cast<double>(report.selected_aggregate_h1_norm2)
        << ", \"normalization_common_factor\": "
        << static_cast<double>(report.normalization_common_factor)
        << ", \"squared_cauchy_majorant\": "
        << static_cast<double>(report.squared_cauchy_majorant)
        << ", \"diagonal_squared_majorant\": "
        << static_cast<double>(report.diagonal_squared_majorant)
        << ", \"schur_squared_majorant\": "
        << static_cast<double>(report.schur_squared_majorant)
        << ", \"height_half_compensated_squared_cauchy_majorant\": "
        << static_cast<double>(
               report.height_half_compensated_squared_cauchy_majorant)
        << ", \"height_half_compensated_diagonal_squared_majorant\": "
        << static_cast<double>(
               report.height_half_compensated_diagonal_squared_majorant)
        << ", \"height_half_compensated_schur_squared_majorant\": "
        << static_cast<double>(
               report.height_half_compensated_schur_squared_majorant)
        << ", \"height_half_compensated_half_decay_squared_majorant\": "
        << static_cast<double>(
               report.height_half_compensated_half_decay_squared_majorant)
        << ", \"exact_tail_h2_reconstruction\": "
        << (report.exact_tail_h2_reconstruction ? "true" : "false")
        << ", \"finite_schur_inequality_verified\": "
        << (report.finite_schur_inequality_verified ? "true" : "false")
        << ", \"finite_half_gap_decay_inequality_verified\": "
        << (report.finite_half_gap_decay_inequality_verified
                ? "true" : "false")
        << ", \"finite\": " << (report.finite ? "true" : "false")
        << ", \"shells\": [\n";
    for (std::size_t index = 0; index < report.shells.size(); ++index) {
        const auto& shell = report.shells[index];
        output << "      {\"shell\": " << shell.shell
            << ", \"minimum_primitive_height\": "
            << shell.minimum_primitive_height
            << ", \"maximum_primitive_height\": "
            << shell.maximum_primitive_height
            << ", \"shape_count\": " << shell.shape_count
            << ", \"interaction_count\": " << shell.interaction_count
            << ", \"aggregate_h2_norm2\": "
            << static_cast<double>(shell.aggregate_h2_norm2)
            << ", \"diagonal_tail_fraction\": "
            << static_cast<double>(shell.diagonal_tail_fraction)
            << ", \"normalized_absolute_gram_row_sum\": "
            << static_cast<double>(
                   shell.normalized_absolute_gram_row_sum)
            << '}'
            << (index + 1 == report.shells.size() ? "\n" : ",\n");
    }
    output << "    ], \"pairs\": [\n";
    for (std::size_t index = 0; index < report.pairs.size(); ++index) {
        const auto& pair = report.pairs[index];
        output << "      {\"first_shell\": " << pair.first_shell
            << ", \"second_shell\": " << pair.second_shell
            << ", \"shell_gap\": " << pair.shell_gap
            << ", \"gram_pairing\": "
            << static_cast<double>(pair.gram_pairing)
            << ", \"normalized_absolute_gram\": "
            << static_cast<double>(pair.normalized_absolute_gram)
            << ", \"half_decay_weighted_correlation\": "
            << static_cast<double>(
                   pair.half_decay_weighted_correlation)
            << '}'
            << (index + 1 == report.pairs.size() ? "\n" : ",\n");
    }
    output << "    ], \"gaps\": [\n";
    for (std::size_t index = 0; index < report.gaps.size(); ++index) {
        const auto& gap = report.gaps[index];
        output << "      {\"shell_gap\": " << gap.shell_gap
            << ", \"pair_count\": " << gap.pair_count
            << ", \"maximum_normalized_absolute_gram\": "
            << static_cast<double>(
                   gap.maximum_normalized_absolute_gram)
            << ", \"sum_normalized_absolute_gram\": "
            << static_cast<double>(gap.sum_normalized_absolute_gram)
            << ", \"maximum_half_decay_weighted_correlation\": "
            << static_cast<double>(
                   gap.maximum_half_decay_weighted_correlation)
            << '}'
            << (index + 1 == report.gaps.size() ? "\n" : ",\n");
    }
    output << "    ]}";
}

}  // namespace

LocalSldProjectiveNormalizationTailSchurCliOptions
LocalSldProjectiveNormalizationTailSchurCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationTailSchurCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--height") {
            options.core_maximum_heights.push_back(
                static_cast<SpectralInteger>(
                    std::stoll(next(index, name))));
        } else if (name == "--state") {
            options.state_paths.push_back(next(index, name));
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown normalization-tail-Schur option: " + name);
        }
    }
    if (options.core_maximum_heights.empty() ||
        options.core_maximum_heights.size() != options.state_paths.size() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.threads < 1 || options.threads > 256 ||
        std::any_of(
            options.core_maximum_heights.begin(),
            options.core_maximum_heights.end(),
            [](SpectralInteger height) { return height < 1; })) {
        throw std::invalid_argument(
            "normalization-tail-Schur requires paired positive heights/states, a certificate, valid selection, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveNormalizationTailSchurCli::print_help(
    std::ostream& out) {
    out << "Projective normalization tail Gram/Schur options:\n"
        << "  --height H          add a core height; repeatable and paired by order\n"
        << "  --state PATH        add its state; repeatable and paired by order\n"
        << "  --selection NAME    local SLD triad selection\n"
        << "  --threads N         parallel direct workers\n"
        << "  --certificate PATH  write English JSON ledger\n";
}

int LocalSldProjectiveNormalizationTailSchurCli::run(
    const LocalSldProjectiveNormalizationTailSchurCliOptions& options,
    std::ostream& out) {
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    std::vector<LocalSldProjectiveNormalizationTailSchurReport> reports;
    reports.reserve(options.state_paths.size());
    bool every_report_finite = true;
    for (std::size_t index = 0; index < options.state_paths.size(); ++index) {
        const SpectralState state = SpectralStateReader::read_tsv(
            options.state_paths[index]);
        reports.push_back(
            LocalSldProjectiveNormalizationTailSchurLedger::analyze(
                dynamics, state,
                LocalSldTriadSelection::parse(options.selection),
                options.core_maximum_heights[index], options.threads));
        every_report_finite =
            every_report_finite && reports.back().finite;
    }
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write normalization-tail-Schur certificate");
    }
    certificate << std::setprecision(18);
    certificate
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-tail-schur-v1\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"exact_identity\": \"T2=sum_(i,j>H)<A b_i,A b_j>\",\n"
        << "  \"finite_schur_bound\": \"T2 <= R_H sum_(i>H)||A b_i||_2^2\",\n"
        << "  \"fitted_squared_cauchy_height_slope\": "
        << static_cast<double>(fitted_height_slope(
               reports, &Report::squared_cauchy_majorant)) << ",\n"
        << "  \"fitted_diagonal_squared_height_slope\": "
        << static_cast<double>(fitted_height_slope(
               reports, &Report::diagonal_squared_majorant)) << ",\n"
        << "  \"fitted_schur_squared_height_slope\": "
        << static_cast<double>(fitted_height_slope(
               reports, &Report::schur_squared_majorant)) << ",\n"
        << "  \"maximum_normalized_absolute_gram_row_sum\": "
        << static_cast<double>(maximum_value(
               reports,
               &Report::maximum_normalized_absolute_gram_row_sum))
        << ",\n"
        << "  \"maximum_half_decay_weighted_correlation\": "
        << static_cast<double>(maximum_value(
               reports,
               &Report::maximum_half_decay_weighted_correlation))
        << ",\n"
        << "  \"maximum_height_half_compensated_diagonal_squared_majorant\": "
        << static_cast<double>(maximum_value(
               reports,
               &Report::height_half_compensated_diagonal_squared_majorant))
        << ",\n"
        << "  \"maximum_height_half_compensated_schur_squared_majorant\": "
        << static_cast<double>(maximum_value(
               reports,
               &Report::height_half_compensated_schur_squared_majorant))
        << ",\n"
        << "  \"maximum_height_half_compensated_half_decay_squared_majorant\": "
        << static_cast<double>(maximum_value(
               reports,
               &Report::height_half_compensated_half_decay_squared_majorant))
        << ",\n"
        << "  \"reports\": [\n";
    for (std::size_t index = 0; index < reports.size(); ++index) {
        write_report(reports[index], options.state_paths[index], certificate);
        certificate << (index + 1 == reports.size() ? "\n" : ",\n");
    }
    certificate
        << "  ],\n"
        << "  \"every_report_finite\": "
        << (every_report_finite ? "true" : "false") << ",\n"
        << "  \"finite_ledger_is_not_a_proof\": true,\n"
        << "  \"uniform_schur_row_bound_proved\": false,\n"
        << "  \"uniform_half_gap_decay_proved\": false,\n"
        << "  \"uniform_diagonal_height_decay_proved\": false,\n"
        << "  \"candidate_split\": \"a joint uniform bound on H^(1/2) times the common normalization times R_H times diagonal tail H2 mass is sufficient for PNT-8\"\n"
        << "}\n";
    out << "projective normalization tail Schur reports="
        << reports.size()
        << " finite=" << (every_report_finite ? "true" : "false")
        << "\nCertificate written to " << options.certificate_path << '\n';
    return every_report_finite ? 0 : 2;
}

}  // namespace lemma
