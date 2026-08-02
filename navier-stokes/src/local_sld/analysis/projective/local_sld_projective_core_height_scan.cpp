#include "local_sld_projective_core_height_scan.hpp"

#include "projective_core_family.hpp"
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

std::vector<SpectralInteger> parse_heights(const std::string& text) {
    std::vector<SpectralInteger> result;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            throw std::invalid_argument(
                "projective core heights contain an empty entry");
        }
        result.push_back(static_cast<SpectralInteger>(std::stoll(token)));
    }
    if (result.empty()) {
        throw std::invalid_argument(
            "projective core heights must be nonempty");
    }
    return result;
}

SpectralReal fit_slope(
    const std::vector<LocalSldProjectiveCoreHeightScanRow>& rows,
    bool absolute_sum) {
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    std::size_t count = 0;
    for (const auto& row : rows) {
        const SpectralReal value = absolute_sum
            ? row.open_absolute_power_one
            : std::abs(row.open_signed_power_one);
        if (row.maximum_height < 1 || value <= 1e-30L) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.maximum_height));
        const SpectralReal y = std::log(value);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        ++count;
    }
    if (count < 2) {
        return 0.0L;
    }
    const SpectralReal n = static_cast<SpectralReal>(count);
    const SpectralReal denominator = n * sum_xx - sum_x * sum_x;
    if (std::abs(denominator) <= 1e-30L) {
        return 0.0L;
    }
    return (n * sum_xy - sum_x * sum_y) / denominator;
}

void write_json(
    const LocalSldProjectiveCoreHeightScanReport& report,
    const LocalSldProjectiveCoreHeightScanCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-core-height-scan-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"fitted_open_signed_height_slope\": "
        << static_cast<double>(report.fitted_open_signed_height_slope)
        << ",\n"
        << "  \"fitted_open_absolute_height_slope\": "
        << static_cast<double>(report.fitted_open_absolute_height_slope)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"maximum_height\": " << row.maximum_height
            << ", \"fixed_family_cardinality\": "
            << row.fixed_family_cardinality
            << ", \"active_core_shape_count\": "
            << row.active_core_shape_count
            << ", \"tail_shape_count\": " << row.tail_shape_count
            << ", \"selected_power_one\": "
            << static_cast<double>(row.selected_power_one)
            << ", \"core_power_one\": "
            << static_cast<double>(row.core_power_one)
            << ", \"core_tail_power_one\": "
            << static_cast<double>(row.core_tail_power_one)
            << ", \"tail_power_one\": "
            << static_cast<double>(row.tail_power_one)
            << ", \"open_signed_power_one\": "
            << static_cast<double>(row.open_signed_power_one)
            << ", \"open_absolute_power_one\": "
            << static_cast<double>(row.open_absolute_power_one)
            << ", \"open_signed_fraction\": "
            << static_cast<double>(row.open_signed_fraction)
            << ", \"open_absolute_fraction\": "
            << static_cast<double>(row.open_absolute_fraction)
            << ", \"decomposition_error\": "
            << static_cast<double>(row.decomposition_error) << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_decomposition_exact\": "
        << (report.every_decomposition_exact ? "true" : "false")
        << ",\n"
        << "  \"each_core_is_fixed_and_finite\": true,\n"
        << "  \"uniform_height_tail_bound_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"interpretation\": \"height slopes are finite-state diagnostics; a summable analytic decay in primitive shape height is still required\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveCoreHeightScanReport
LocalSldProjectiveCoreHeightScan::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    std::vector<SpectralInteger> maximum_heights,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (maximum_heights.empty()) {
        throw std::invalid_argument(
            "projective core-height scan requires heights");
    }
    std::sort(maximum_heights.begin(), maximum_heights.end());
    if (std::adjacent_find(
            maximum_heights.begin(), maximum_heights.end()) !=
        maximum_heights.end()) {
        throw std::invalid_argument(
            "projective core-height scan contains duplicate heights");
    }
    LocalSldProjectiveCoreHeightScanReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.every_decomposition_exact = true;
    report.each_core_is_fixed_and_finite = true;
    for (const SpectralInteger height : maximum_heights) {
        const std::vector<ProjectivePrimitiveSignature> core =
            ProjectiveCoreFamily::through_maximum_height(height);
        const LocalSldProjectiveCoreTailReport decomposition =
            LocalSldProjectiveCoreTailLedger::analyze(
                dynamics, state, core, threads,
                exclude_signature_123, exclude_triple_family);
        LocalSldProjectiveCoreHeightScanRow row;
        row.maximum_height = height;
        row.fixed_family_cardinality = core.size();
        row.active_core_shape_count =
            decomposition.active_core_shape_count;
        row.tail_shape_count = decomposition.tail_shape_count;
        row.selected_power_one = decomposition.selected_power_one;
        row.core_power_one = decomposition.core.power_one;
        row.core_tail_power_one = decomposition.core_tail.power_one;
        row.tail_power_one = decomposition.tail.power_one;
        row.open_signed_power_one =
            row.core_tail_power_one + row.tail_power_one;
        row.open_absolute_power_one =
            std::abs(row.core_tail_power_one) +
            std::abs(row.tail_power_one);
        const SpectralReal denominator = std::max(
            std::abs(row.selected_power_one), 1e-30L);
        row.open_signed_fraction =
            std::abs(row.open_signed_power_one) / denominator;
        row.open_absolute_fraction =
            row.open_absolute_power_one / denominator;
        row.decomposition_error = decomposition.bracket_partition_error;
        report.every_decomposition_exact =
            report.every_decomposition_exact &&
            decomposition.exact_core_tail_decomposition;
        report.rows.push_back(row);
    }
    report.fitted_open_signed_height_slope = fit_slope(
        report.rows, false);
    report.fitted_open_absolute_height_slope = fit_slope(
        report.rows, true);
    return report;
}

LocalSldProjectiveCoreHeightScanCliOptions
LocalSldProjectiveCoreHeightScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveCoreHeightScanCliOptions options;
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
        } else if (name == "--heights") {
            options.maximum_heights = parse_heights(next(index, name));
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
                "unknown projective-core-height-scan option: " + name);
        }
    }
    if (options.state_path.empty() || options.maximum_heights.empty() ||
        options.certificate_path.empty() || options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "projective-core-height-scan requires state, heights, certificate, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveCoreHeightScanCli::print_help(
    std::ostream& out) {
    out << "Projective primitive-shape-height core scan options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --heights A,B,C       increasing fixed core height thresholds\n"
        << "  --certificate PATH    write English JSON scan\n"
        << "  --threads N           parallel projective interaction workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveCoreHeightScanCli::run(
    const LocalSldProjectiveCoreHeightScanCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldProjectiveCoreHeightScanReport report =
        LocalSldProjectiveCoreHeightScan::analyze(
            dynamics, state, options.maximum_heights,
            options.threads, options.exclude_signature_123,
            options.exclude_triple_family);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective core-height scan certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective core-height scan cutoff=" << report.cutoff
        << " rows=" << report.rows.size()
        << " signed_height_slope="
        << static_cast<double>(report.fitted_open_signed_height_slope)
        << " absolute_height_slope="
        << static_cast<double>(report.fitted_open_absolute_height_slope)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_decomposition_exact ? 0 : 2;
}

}  // namespace lemma
