#include "local_sld_projective_height_coercivity_line_scan.hpp"

#include "local_sld_projective_height_commutator_ratio_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_state_blend.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

SpectralReal low_outer_exponent(
    const std::vector<LocalSldProjectiveHeightCoercivityLineRow>& rows) {
    std::vector<std::size_t> order(rows.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t left,
                                               std::size_t right) {
        return rows[left].outer_h1_weight < rows[right].outer_h1_weight;
    });
    const std::size_t count = std::min<std::size_t>(7, order.size());
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    std::size_t valid = 0;
    for (std::size_t position = 0; position < count; ++position) {
        const auto& row = rows[order[position]];
        if (!(row.outer_h1_weight > 0.0L) ||
            !(row.paired_envelope > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(row.outer_h1_weight);
        const SpectralReal y = std::log(row.paired_envelope);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        ++valid;
    }
    if (valid < 2) {
        return std::numeric_limits<SpectralReal>::quiet_NaN();
    }
    const SpectralReal n = static_cast<SpectralReal>(valid);
    const SpectralReal denominator = n * sum_xx - sum_x * sum_x;
    return std::abs(denominator) > 1e-30L
        ? (n * sum_xy - sum_x * sum_y) / denominator
        : std::numeric_limits<SpectralReal>::quiet_NaN();
}

void write_json(
    const LocalSldProjectiveHeightCoercivityLineReport& report,
    const LocalSldProjectiveHeightCoercivityLineScanOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-coercivity-line-v1\",\n"
        << "  \"left_state_path\": \"" << options.left_state_path
        << "\",\n"
        << "  \"right_state_path\": \"" << options.right_state_path
        << "\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"path_formula\": \"normalize((1-t) left + t right)\",\n"
        << "  \"maximum_ratio\": "
        << static_cast<double>(report.maximum_ratio) << ",\n"
        << "  \"maximum_ratio_parameter\": "
        << static_cast<double>(report.maximum_ratio_parameter) << ",\n"
        << "  \"minimum_outer_h1_weight\": "
        << static_cast<double>(report.minimum_outer_h1_weight) << ",\n"
        << "  \"minimum_outer_parameter\": "
        << static_cast<double>(report.minimum_outer_parameter) << ",\n"
        << "  \"low_outer_paired_exponent\": "
        << static_cast<double>(report.low_outer_paired_exponent) << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"parameter\": "
            << static_cast<double>(row.parameter)
            << ", \"paired_envelope\": "
            << static_cast<double>(row.paired_envelope)
            << ", \"outer_h1_weight\": "
            << static_cast<double>(row.outer_h1_weight)
            << ", \"coercivity_ratio\": "
            << static_cast<double>(row.coercivity_ratio)
            << ", \"squared_coercivity_ratio\": "
            << static_cast<double>(row.squared_coercivity_ratio)
            << "}" << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"finite_line_scan_is_not_a_proof\": true,\n"
        << "  \"candidate_outer_coercivity_lemma_proved\": false\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveHeightCoercivityLineReport
LocalSldProjectiveHeightCoercivityLineScan::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& left,
    const SpectralState& right,
    TriadSelection selection,
    SpectralReal minimum_parameter,
    SpectralReal maximum_parameter,
    int samples,
    int threads) {
    if (!std::isfinite(minimum_parameter) ||
        !std::isfinite(maximum_parameter) ||
        !(maximum_parameter > minimum_parameter) || samples < 2 ||
        samples > 1000 || threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "invalid projective-height coercivity line parameters");
    }
    const LocalSldProjectiveHeightCommutatorRatioObjective objective(
        dynamics, selection, threads);
    LocalSldProjectiveHeightCoercivityLineReport report;
    report.cutoff = std::max(
        SpectralStateOps::cutoff(left), SpectralStateOps::cutoff(right));
    report.minimum_outer_h1_weight =
        std::numeric_limits<SpectralReal>::infinity();
    report.rows.reserve(static_cast<std::size_t>(samples));
    for (int index = 0; index < samples; ++index) {
        const SpectralReal fraction = static_cast<SpectralReal>(index) /
            static_cast<SpectralReal>(samples - 1);
        const SpectralReal parameter = minimum_parameter +
            fraction * (maximum_parameter - minimum_parameter);
        const SpectralState state = SpectralStateBlend::affine_normalized(
            left, right, parameter);
        const auto value = objective.evaluate(state);
        LocalSldProjectiveHeightCoercivityLineRow row;
        row.parameter = parameter;
        row.paired_envelope = value.commutator_paired_bracket_envelope;
        row.outer_h1_weight = value.outer_h1_sum;
        row.coercivity_ratio = value.coercivity_ratio;
        row.squared_coercivity_ratio = value.squared_coercivity_ratio;
        if (row.coercivity_ratio > report.maximum_ratio) {
            report.maximum_ratio = row.coercivity_ratio;
            report.maximum_ratio_parameter = parameter;
        }
        if (row.outer_h1_weight < report.minimum_outer_h1_weight) {
            report.minimum_outer_h1_weight = row.outer_h1_weight;
            report.minimum_outer_parameter = parameter;
        }
        report.rows.push_back(row);
    }
    report.low_outer_paired_exponent = low_outer_exponent(report.rows);
    return report;
}

LocalSldProjectiveHeightCoercivityLineScanOptions
LocalSldProjectiveHeightCoercivityLineScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightCoercivityLineScanOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--left-state") {
            options.left_state_path = next(index, name);
        } else if (name == "--right-state") {
            options.right_state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--selection") {
            options.selection = next(index, name);
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
                "unknown projective-height coercivity-line option: " +
                name);
        }
    }
    if (options.left_state_path.empty() ||
        options.right_state_path.empty() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        !std::isfinite(options.minimum_parameter) ||
        !std::isfinite(options.maximum_parameter) ||
        !(options.maximum_parameter > options.minimum_parameter) ||
        options.samples < 2 || options.samples > 1000 ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "coercivity-line requires two states, certificate, valid selection and parameter range, samples 2..1000, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveHeightCoercivityLineScanCli::print_help(
    std::ostream& out) {
    out << "Projective-height outer-coercivity affine-line options:\n"
        << "  --left-state PATH     affine state at t=0\n"
        << "  --right-state PATH    affine state at t=1\n"
        << "  --certificate PATH    write English JSON scan\n"
        << "  --selection NAME      local SLD triad selection\n"
        << "  --min-parameter X     first affine parameter\n"
        << "  --max-parameter X     last affine parameter\n"
        << "  --samples N           equally spaced samples\n"
        << "  --threads N           parallel projective workers\n";
}

int LocalSldProjectiveHeightCoercivityLineScanCli::run(
    const LocalSldProjectiveHeightCoercivityLineScanOptions& options,
    std::ostream& out) {
    const SpectralState left = SpectralStateReader::read_tsv(
        options.left_state_path);
    const SpectralState right = SpectralStateReader::read_tsv(
        options.right_state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const auto report = LocalSldProjectiveHeightCoercivityLineScan::analyze(
        dynamics, left, right,
        LocalSldTriadSelection::parse(options.selection),
        options.minimum_parameter, options.maximum_parameter,
        options.samples, options.threads);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective-height coercivity-line certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective-height coercivity line cutoff=" << report.cutoff
        << " maximum_ratio="
        << static_cast<double>(report.maximum_ratio)
        << " at=" << static_cast<double>(report.maximum_ratio_parameter)
        << " minimum_outer="
        << static_cast<double>(report.minimum_outer_h1_weight)
        << " at=" << static_cast<double>(report.minimum_outer_parameter)
        << " paired_exponent="
        << static_cast<double>(report.low_outer_paired_exponent)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
