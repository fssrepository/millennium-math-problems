#include "local_sld_projective_fan_scan.hpp"

#include "local_sld_projective_fan_state.hpp"
#include "local_sld_shape_power_objective.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

template <class Getter>
SpectralReal fitted_slope(
    const std::vector<LocalSldProjectiveFanScanRow>& rows,
    Getter getter) {
    SpectralReal n = 0.0L;
    SpectralReal sx = 0.0L;
    SpectralReal sy = 0.0L;
    SpectralReal sxx = 0.0L;
    SpectralReal sxy = 0.0L;
    for (const auto& row : rows) {
        const SpectralReal value = std::abs(getter(row));
        if (row.cutoff < 1 || !(value > 1e-30L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
        const SpectralReal y = std::log(value);
        n += 1.0L;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    const SpectralReal denominator = n * sxx - sx * sx;
    return n >= 2.0L && std::abs(denominator) > 1e-30L
        ? (n * sxy - sx * sy) / denominator
        : 0.0L;
}

}  // namespace

LocalSldProjectiveFanScanReport
LocalSldProjectiveFanScan::analyze(
    const LocalSldProjectiveFanScanOptions& options) {
    if (options.minimum_cutoff < 3 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 24 || options.threads < 1) {
        throw std::invalid_argument(
            "invalid projective fan scan options");
    }
    std::filesystem::create_directories(options.state_directory);
    LocalSldProjectiveFanScanReport report;
    report.threads = options.threads;
    const TriadSelection selection =
        TriadSelection::local_without_equal_low_double_triple();
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        LocalSldProjectiveFanState generated =
            LocalSldProjectiveFanStateFactory::make(cutoff);
        SpectralGalerkin configuration;
        configuration.configure("direct", options.threads);
        const SpectralDynamics dynamics(configuration);
        const LocalSldProjectiveCoherenceReport coherence =
            LocalSldProjectiveCoherenceLedger::analyze(
                dynamics, generated.state, options.threads,
                false, true);
        const LocalSldShapePowerObjectiveValue power_one =
            LocalSldShapePowerObjective(
                dynamics, selection, 1).evaluate(generated.state);
        LocalSldProjectiveFanScanRow row;
        row.cutoff = cutoff;
        row.coherent_pairs = generated.coherent_pairs;
        row.active_positive_modes = generated.active_positive_modes;
        row.projective_shapes = coherence.projective_shape_count;
        row.synthesis_ratio = coherence.coherent_synthesis_ratio;
        row.synthesis_amplification =
            coherence.coherent_synthesis_amplification;
        row.maximum_output_synthesis_ratio =
            coherence.maximum_output_synthesis_ratio;
        row.maximum_output_wave = coherence.maximum_output_wave;
        row.bracket_constant_ratio =
            power_one.bracket_constant_ratio;
        row.normalized_full_stretching =
            power_one.normalized_stretching;
        row.power_one_product = power_one.absolute_power_product;
        row.projective_reconstruction_error =
            coherence.reconstruction_relative_error;
        row.state_path = (
            std::filesystem::path(options.state_directory) /
            ("K" + std::to_string(cutoff) + ".tsv")).string();
        std::ostringstream metadata;
        metadata << "explicit planar projective coherent-fan state; "
            << "cutoff=" << cutoff
            << "; coherent_pairs=" << row.coherent_pairs
            << "; candidate_lemma_proved=false";
        SpectralStateWriter::write_tsv(
            row.state_path, generated.state, metadata.str());
        report.rows.push_back(row);
    }
    report.synthesis_ratio_cutoff_slope = fitted_slope(
        report.rows, [](const auto& row) {
            return row.synthesis_ratio;
        });
    report.power_one_cutoff_slope = fitted_slope(
        report.rows, [](const auto& row) {
            return row.power_one_product;
        });
    report.every_projective_reconstruction_exact = std::all_of(
        report.rows.begin(), report.rows.end(),
        [](const auto& row) {
            return row.projective_reconstruction_error < 1e-13L;
        });
    return report;
}

LocalSldProjectiveFanScanOptions LocalSldProjectiveFanScan::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveFanScanOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--min-cutoff") {
            options.minimum_cutoff = std::stoi(next(index, name));
        } else if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--state-dir") {
            options.state_directory = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown projective-fan-scan option: " + name);
        }
    }
    if (options.minimum_cutoff < 3 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 24 || options.threads < 1 ||
        options.certificate_path.empty() ||
        options.state_directory.empty()) {
        throw std::invalid_argument(
            "projective-fan-scan requires cutoffs 3..24, positive threads, certificate, and state directory");
    }
    return options;
}

void LocalSldProjectiveFanScan::print_help(std::ostream& out) {
    out << "Explicit projective coherent-fan scan options:\n"
        << "  --min-cutoff K       first planar cutoff (minimum 3)\n"
        << "  --max-cutoff K       last planar cutoff (maximum 24)\n"
        << "  --threads N          projective-ray workers\n"
        << "  --certificate PATH   write English JSON scan\n"
        << "  --state-dir PATH     write replayable fan states\n";
}

int LocalSldProjectiveFanScan::run(
    const LocalSldProjectiveFanScanOptions& options,
    std::ostream& out) {
    const LocalSldProjectiveFanScanReport report = analyze(options);
    const std::filesystem::path certificate_path(
        options.certificate_path);
    if (!certificate_path.parent_path().empty()) {
        std::filesystem::create_directories(
            certificate_path.parent_path());
    }
    std::ofstream certificate(certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective fan certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-fan-scan-v1\",\n"
        << "  \"construction\": \"planar x-invariant modes p=(0,-a,K-b), q=(0,a,b), p+q=(0,0,K), with phases chosen so the target contributions align\",\n"
        << "  \"selection\": \"local signatures excluding (m,m,2m) and (m,m,3m)\",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"synthesis_ratio_cutoff_slope\": "
        << static_cast<double>(report.synthesis_ratio_cutoff_slope)
        << ",\n"
        << "  \"power_one_cutoff_slope\": "
        << static_cast<double>(report.power_one_cutoff_slope) << ",\n"
        << "  \"every_projective_reconstruction_exact\": "
        << (report.every_projective_reconstruction_exact
            ? "true" : "false") << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"coherent_pairs\": " << row.coherent_pairs
            << ", \"active_positive_modes\": "
            << row.active_positive_modes
            << ", \"projective_shapes\": " << row.projective_shapes
            << ", \"synthesis_ratio\": "
            << static_cast<double>(row.synthesis_ratio)
            << ", \"synthesis_amplification\": "
            << static_cast<double>(row.synthesis_amplification)
            << ", \"maximum_output_synthesis_ratio\": "
            << static_cast<double>(row.maximum_output_synthesis_ratio)
            << ", \"maximum_output_wave\": ["
            << row.maximum_output_wave.x << ", "
            << row.maximum_output_wave.y << ", "
            << row.maximum_output_wave.z << "]"
            << ", \"bracket_constant_ratio\": "
            << static_cast<double>(row.bracket_constant_ratio)
            << ", \"normalized_full_stretching\": "
            << static_cast<double>(row.normalized_full_stretching)
            << ", \"power_one_product\": "
            << static_cast<double>(row.power_one_product)
            << ", \"projective_reconstruction_error\": "
            << static_cast<double>(
                   row.projective_reconstruction_error)
            << ", \"state_path\": \"" << row.state_path << "\"}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    certificate
        << "  ],\n"
        << "  \"cutoff_independent_synthesis_bound_proved\": false,\n"
        << "  \"power_one_tradeoff_bound_proved\": false,\n"
        << "  \"finite_scan_is_not_a_proof\": true\n"
        << "}\n";
    out << std::setprecision(12);
    for (const auto& row : report.rows) {
        out << "projective fan cutoff=" << row.cutoff
            << " pairs=" << row.coherent_pairs
            << " shapes=" << row.projective_shapes
            << " synthesis=" << static_cast<double>(row.synthesis_ratio)
            << " max_output="
            << static_cast<double>(row.maximum_output_synthesis_ratio)
            << " power_one="
            << static_cast<double>(row.power_one_product) << '\n';
    }
    out << "synthesis slope="
        << static_cast<double>(report.synthesis_ratio_cutoff_slope)
        << " power-one slope="
        << static_cast<double>(report.power_one_cutoff_slope) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_projective_reconstruction_exact ? 0 : 2;
}

}  // namespace lemma
