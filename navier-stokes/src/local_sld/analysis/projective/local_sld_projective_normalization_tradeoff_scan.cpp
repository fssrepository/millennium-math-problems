#include "local_sld_projective_normalization_tradeoff_scan.hpp"

#include "local_sld_projective_core_tail_alignment_objective.hpp"
#include "local_sld_projective_normalization_alignment_objective.hpp"
#include "local_sld_projective_normalization_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_state_blend.hpp"
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

void write_json(
    const LocalSldProjectiveNormalizationTradeoffReport& report,
    const LocalSldProjectiveNormalizationTradeoffScanOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-tradeoff-line-v2\",\n"
        << "  \"normalization_state_path\": \""
        << options.normalization_state_path << "\",\n"
        << "  \"alignment_state_path\": \""
        << options.alignment_state_path << "\",\n"
        << "  \"triad_selection\": \"" << options.selection
        << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"core_maximum_height\": "
        << report.core_maximum_height << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"path_formula\": \"normalize((1-t) normalization_state + t alignment_state)\",\n"
        << "  \"maximum_open_power_one\": "
        << static_cast<double>(report.maximum_open_power_one) << ",\n"
        << "  \"maximum_open_parameter\": "
        << static_cast<double>(report.maximum_open_parameter) << ",\n"
        << "  \"maximum_tail_alignment_squared\": "
        << static_cast<double>(report.maximum_tail_alignment_squared)
        << ",\n"
        << "  \"maximum_tail_alignment_parameter\": "
        << static_cast<double>(report.maximum_tail_alignment_parameter)
        << ",\n"
        << "  \"maximum_normalization_alignment_product_squared\": "
        << static_cast<double>(
               report.maximum_normalization_alignment_product_squared)
        << ",\n"
        << "  \"maximum_normalization_alignment_parameter\": "
        << static_cast<double>(
               report.maximum_normalization_alignment_parameter)
        << ",\n"
        << "  \"maximum_selected_channel_power_one\": "
        << static_cast<double>(report.maximum_selected_channel_power_one)
        << ",\n"
        << "  \"maximum_selected_channel_parameter\": "
        << static_cast<double>(report.maximum_selected_channel_parameter)
        << ",\n"
        << "  \"maximum_roughness\": "
        << static_cast<double>(report.maximum_roughness) << ",\n"
        << "  \"maximum_roughness_parameter\": "
        << static_cast<double>(report.maximum_roughness_parameter)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"parameter\": "
            << static_cast<double>(row.parameter)
            << ", \"open_palinstrophy_normalization_power_one\": "
            << static_cast<double>(
                   row.open_palinstrophy_normalization_power_one)
            << ", \"squared_open_palinstrophy_normalization_power_one\": "
            << static_cast<double>(
                   row.squared_open_palinstrophy_normalization_power_one)
            << ", \"tail_stretching_alignment_squared\": "
            << static_cast<double>(
                   row.tail_stretching_alignment_squared)
            << ", \"selected_stretching_h1_alignment_squared\": "
            << static_cast<double>(
                   row.selected_stretching_h1_alignment_squared)
            << ", \"tail_palinstrophy_cross_h2_alignment_squared\": "
            << static_cast<double>(
                   row.tail_palinstrophy_cross_h2_alignment_squared)
            << ", \"normalization_alignment_product_squared\": "
            << static_cast<double>(
                   row.normalization_alignment_product_squared)
            << ", \"selected_stretching_tail_cross_power_one\": "
            << static_cast<double>(
                   row.selected_stretching_tail_cross_power_one)
            << ", \"enstrophy\": "
            << static_cast<double>(row.enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(row.palinstrophy)
            << ", \"palinstrophy_over_enstrophy_squared\": "
            << static_cast<double>(
                   row.palinstrophy_over_enstrophy_squared)
            << ", \"full_stretching\": "
            << static_cast<double>(row.full_stretching)
            << ", \"normalization_scale\": "
            << static_cast<double>(row.normalization_scale)
            << ", \"unscaled_open_normalization\": "
            << static_cast<double>(row.unscaled_open_normalization)
            << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_sample_finite\": "
        << (report.every_sample_finite ? "true" : "false") << ",\n"
        << "  \"finite_line_scan_is_not_a_proof\": true,\n"
        << "  \"uniform_scale_alignment_tradeoff_proved\": false,\n"
        << "  \"interpretation\": \"large tail alignment can coexist with large unscaled interaction, but the Z^2 P^2 normalization suppresses the rough endpoint; the missing lemma must control their joint scale-aware product\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveNormalizationTradeoffReport
LocalSldProjectiveNormalizationTradeoffScan::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& normalization_state,
    const SpectralState& alignment_state,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    SpectralReal minimum_parameter,
    SpectralReal maximum_parameter,
    int samples,
    int threads) {
    if (core_maximum_height < 1 ||
        !std::isfinite(minimum_parameter) ||
        !std::isfinite(maximum_parameter) ||
        !(maximum_parameter > minimum_parameter) ||
        samples < 2 || samples > 1000 ||
        threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "invalid normalization tradeoff-line parameters");
    }
    const LocalSldProjectiveNormalizationObjective normalization(
        dynamics, selection, core_maximum_height, threads);
    const LocalSldProjectiveCoreTailAlignmentObjective alignment(
        dynamics, selection, core_maximum_height,
        LocalSldProjectiveHeightRegion::tail, threads);
    const LocalSldProjectiveNormalizationAlignmentObjective
        normalization_alignment(
            dynamics, selection, core_maximum_height, threads);
    LocalSldProjectiveNormalizationTradeoffReport report;
    report.cutoff = std::max(
        SpectralStateOps::cutoff(normalization_state),
        SpectralStateOps::cutoff(alignment_state));
    report.core_maximum_height = core_maximum_height;
    report.every_sample_finite = true;
    report.rows.reserve(static_cast<std::size_t>(samples));
    for (int index = 0; index < samples; ++index) {
        const SpectralReal fraction = static_cast<SpectralReal>(index) /
            static_cast<SpectralReal>(samples - 1);
        const SpectralReal parameter = minimum_parameter +
            fraction * (maximum_parameter - minimum_parameter);
        const SpectralState state = SpectralStateBlend::affine_normalized(
            normalization_state, alignment_state, parameter);
        const auto normalization_value = normalization.evaluate(state);
        const auto alignment_value = alignment.evaluate(state);
        const auto normalization_alignment_value =
            normalization_alignment.evaluate(state);
        LocalSldProjectiveNormalizationTradeoffRow row;
        row.parameter = parameter;
        row.open_palinstrophy_normalization_power_one =
            normalization_value.palinstrophy_normalization_power_one;
        row.squared_open_palinstrophy_normalization_power_one =
            normalization_value
                .squared_palinstrophy_normalization_power_one;
        row.tail_stretching_alignment_squared =
            alignment_value.stretching_h1_alignment_squared;
        row.selected_stretching_h1_alignment_squared =
            normalization_alignment_value
                .selected_stretching_h1_alignment_squared;
        row.tail_palinstrophy_cross_h2_alignment_squared =
            normalization_alignment_value
                .tail_palinstrophy_cross_h2_alignment_squared;
        row.normalization_alignment_product_squared =
            normalization_alignment_value
                .normalization_alignment_product_squared;
        row.selected_stretching_tail_cross_power_one =
            normalization_value
                .selected_stretching_tail_cross_power_one;
        row.enstrophy = normalization_value.enstrophy;
        row.palinstrophy = normalization_value.palinstrophy;
        if (row.enstrophy > 0.0L) {
            row.palinstrophy_over_enstrophy_squared =
                row.palinstrophy /
                (row.enstrophy * row.enstrophy);
        }
        row.full_stretching = normalization_value.full_stretching;
        if (row.enstrophy > 0.0L && row.palinstrophy > 0.0L) {
            row.normalization_scale =
                std::abs(row.full_stretching) /
                (row.enstrophy * row.enstrophy *
                 row.palinstrophy * row.palinstrophy);
        }
        row.unscaled_open_normalization =
            normalization_value.open_palinstrophy_normalization;
        report.every_sample_finite = report.every_sample_finite &&
            normalization_value.finite && alignment_value.finite &&
            normalization_alignment_value.finite &&
            std::isfinite(row.palinstrophy_over_enstrophy_squared) &&
            std::isfinite(row.normalization_scale);
        if (row.open_palinstrophy_normalization_power_one >
            report.maximum_open_power_one) {
            report.maximum_open_power_one =
                row.open_palinstrophy_normalization_power_one;
            report.maximum_open_parameter = parameter;
        }
        if (row.tail_stretching_alignment_squared >
            report.maximum_tail_alignment_squared) {
            report.maximum_tail_alignment_squared =
                row.tail_stretching_alignment_squared;
            report.maximum_tail_alignment_parameter = parameter;
        }
        if (row.normalization_alignment_product_squared >
            report.maximum_normalization_alignment_product_squared) {
            report.maximum_normalization_alignment_product_squared =
                row.normalization_alignment_product_squared;
            report.maximum_normalization_alignment_parameter = parameter;
        }
        if (row.selected_stretching_tail_cross_power_one >
            report.maximum_selected_channel_power_one) {
            report.maximum_selected_channel_power_one =
                row.selected_stretching_tail_cross_power_one;
            report.maximum_selected_channel_parameter = parameter;
        }
        if (row.palinstrophy_over_enstrophy_squared >
            report.maximum_roughness) {
            report.maximum_roughness =
                row.palinstrophy_over_enstrophy_squared;
            report.maximum_roughness_parameter = parameter;
        }
        report.rows.push_back(row);
    }
    return report;
}

LocalSldProjectiveNormalizationTradeoffScanOptions
LocalSldProjectiveNormalizationTradeoffScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationTradeoffScanOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--normalization-state") {
            options.normalization_state_path = next(index, name);
        } else if (name == "--alignment-state") {
            options.alignment_state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--projective-core-height") {
            options.core_maximum_height =
                static_cast<SpectralInteger>(
                    std::stoll(next(index, name)));
        } else if (name == "--min-parameter") {
            options.minimum_parameter = std::stold(next(index, name));
        } else if (name == "--max-parameter") {
            options.maximum_parameter = std::stold(next(index, name));
        } else if (name == "--samples") {
            options.samples = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown normalization-tradeoff option: " + name);
        }
    }
    if (options.normalization_state_path.empty() ||
        options.alignment_state_path.empty() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.core_maximum_height < 1 ||
        !std::isfinite(options.minimum_parameter) ||
        !std::isfinite(options.maximum_parameter) ||
        !(options.maximum_parameter > options.minimum_parameter) ||
        options.samples < 2 || options.samples > 1000 ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "normalization-tradeoff requires two states, certificate, valid selection/core height, parameter range, samples 2..1000, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveNormalizationTradeoffScanCli::print_help(
    std::ostream& out) {
    out << "Projective normalization/alignment affine tradeoff options:\n"
        << "  --normalization-state PATH  open-normalization endpoint\n"
        << "  --alignment-state PATH      tail-alignment endpoint\n"
        << "  --projective-core-height H  fixed primitive-height core\n"
        << "  --selection NAME            local SLD triad selection\n"
        << "  --min-parameter X           first affine parameter\n"
        << "  --max-parameter X           last affine parameter\n"
        << "  --samples N                 equally spaced samples\n"
        << "  --threads N                 parallel projective workers\n"
        << "  --certificate PATH          write English JSON scan\n";
}

int LocalSldProjectiveNormalizationTradeoffScanCli::run(
    const LocalSldProjectiveNormalizationTradeoffScanOptions& options,
    std::ostream& out) {
    const SpectralState normalization_state =
        SpectralStateReader::read_tsv(options.normalization_state_path);
    const SpectralState alignment_state =
        SpectralStateReader::read_tsv(options.alignment_state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const auto report = LocalSldProjectiveNormalizationTradeoffScan::analyze(
        dynamics, normalization_state, alignment_state,
        LocalSldTriadSelection::parse(options.selection),
        options.core_maximum_height,
        options.minimum_parameter, options.maximum_parameter,
        options.samples, options.threads);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write normalization tradeoff certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective normalization tradeoff cutoff=" << report.cutoff
        << " max_open="
        << static_cast<double>(report.maximum_open_power_one)
        << " at=" << static_cast<double>(report.maximum_open_parameter)
        << " max_alignment_squared="
        << static_cast<double>(report.maximum_tail_alignment_squared)
        << " at="
        << static_cast<double>(report.maximum_tail_alignment_parameter)
        << " max_joint_alignment_squared="
        << static_cast<double>(
               report.maximum_normalization_alignment_product_squared)
        << " at="
        << static_cast<double>(
               report.maximum_normalization_alignment_parameter)
        << " max_selected_channel="
        << static_cast<double>(report.maximum_selected_channel_power_one)
        << " at="
        << static_cast<double>(report.maximum_selected_channel_parameter)
        << " max_roughness="
        << static_cast<double>(report.maximum_roughness)
        << " at="
        << static_cast<double>(report.maximum_roughness_parameter)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_sample_finite ? 0 : 2;
}

}  // namespace lemma
