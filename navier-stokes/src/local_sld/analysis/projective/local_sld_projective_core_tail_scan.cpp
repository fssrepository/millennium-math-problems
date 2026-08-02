#include "local_sld_projective_core_tail_scan.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"
#include "projective_core_family.hpp"

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

using Signature = LocalSldProjectiveCoreSignature;

Signature parse_signature(const std::string& text) {
    Signature result{};
    std::stringstream stream(text);
    std::string token;
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (!std::getline(stream, token, ',')) {
            throw std::invalid_argument(
                "core signature must have three comma-separated integers");
        }
        result[index] = static_cast<SpectralInteger>(std::stoll(token));
    }
    if (std::getline(stream, token, ',')) {
        throw std::invalid_argument(
            "core signature must have exactly three entries");
    }
    std::sort(result.begin(), result.end());
    return result;
}

SpectralReal fit_slope(
    const std::vector<LocalSldProjectiveCoreTailScanRow>& rows,
    bool absolute_sum,
    std::size_t minimum_active_core_shape_count = 0) {
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    std::size_t count = 0;
    for (const auto& row : rows) {
        if (row.active_core_shape_count <
            minimum_active_core_shape_count) {
            continue;
        }
        const SpectralReal value = absolute_sum
            ? row.open_absolute_power_one
            : std::abs(row.open_signed_power_one);
        if (row.cutoff < 1 || value <= 1e-30L) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
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
    const LocalSldProjectiveCoreTailScanReport& report,
    const LocalSldProjectiveCoreTailScanCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-core-tail-scan-v1\",\n"
        << "  \"state_paths\": [\n";
    for (std::size_t index = 0; index < options.state_paths.size();
         ++index) {
        output << "    \"" << options.state_paths[index] << "\""
            << (index + 1 == options.state_paths.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"fitted_open_signed_cutoff_slope\": "
        << static_cast<double>(report.fitted_open_signed_cutoff_slope)
        << ",\n"
        << "  \"fitted_open_absolute_cutoff_slope\": "
        << static_cast<double>(report.fitted_open_absolute_cutoff_slope)
        << ",\n"
        << "  \"full_active_core_start_cutoff\": "
        << report.full_active_core_start_cutoff << ",\n"
        << "  \"fitted_full_active_core_open_signed_cutoff_slope\": "
        << static_cast<double>(
               report.fitted_full_active_core_open_signed_cutoff_slope)
        << ",\n"
        << "  \"fitted_full_active_core_open_absolute_cutoff_slope\": "
        << static_cast<double>(
               report.fitted_full_active_core_open_absolute_cutoff_slope)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"cutoff\": " << row.cutoff
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
        << "  \"fixed_core_internal_bound_proved\": true,\n"
        << "  \"uniform_open_bound_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"interpretation\": \"the fitted slopes describe only the supplied finite states; the open signed core-tail plus tail-tail estimate remains unproved\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveCoreTailScanReport
LocalSldProjectiveCoreTailScan::analyze(
    const SpectralDynamics& dynamics,
    const std::vector<SpectralState>& states,
    const std::vector<LocalSldProjectiveCoreSignature>& core,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (states.empty()) {
        throw std::invalid_argument(
            "projective core-tail scan requires at least one state");
    }
    LocalSldProjectiveCoreTailScanReport report;
    report.every_decomposition_exact = true;
    report.fixed_core_internal_bound_proved = true;
    for (const SpectralState& state : states) {
        const LocalSldProjectiveCoreTailReport decomposition =
            LocalSldProjectiveCoreTailLedger::analyze(
                dynamics, state, core, threads,
                exclude_signature_123, exclude_triple_family);
        LocalSldProjectiveCoreTailScanRow row;
        row.cutoff = decomposition.cutoff;
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
        row.decomposition_error =
            decomposition.bracket_partition_error;
        report.every_decomposition_exact =
            report.every_decomposition_exact &&
            decomposition.exact_core_tail_decomposition;
        report.rows.push_back(row);
    }
    std::sort(report.rows.begin(), report.rows.end(),
        [](const auto& left, const auto& right) {
            return left.cutoff < right.cutoff;
        });
    for (std::size_t index = 1; index < report.rows.size(); ++index) {
        if (report.rows[index - 1].cutoff == report.rows[index].cutoff) {
            throw std::invalid_argument(
                "projective core-tail scan has duplicate cutoffs");
        }
    }
    report.fitted_open_signed_cutoff_slope = fit_slope(
        report.rows, false);
    report.fitted_open_absolute_cutoff_slope = fit_slope(
        report.rows, true);
    std::size_t maximum_active_core_shape_count = 0;
    for (const auto& row : report.rows) {
        maximum_active_core_shape_count = std::max(
            maximum_active_core_shape_count,
            row.active_core_shape_count);
    }
    for (const auto& row : report.rows) {
        if (row.active_core_shape_count ==
            maximum_active_core_shape_count) {
            report.full_active_core_start_cutoff = row.cutoff;
            break;
        }
    }
    report.fitted_full_active_core_open_signed_cutoff_slope = fit_slope(
        report.rows, false, maximum_active_core_shape_count);
    report.fitted_full_active_core_open_absolute_cutoff_slope = fit_slope(
        report.rows, true, maximum_active_core_shape_count);
    return report;
}

LocalSldProjectiveCoreTailScanCliOptions
LocalSldProjectiveCoreTailScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveCoreTailScanCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_paths.push_back(next(index, name));
        } else if (name == "--core-signature") {
            options.core.push_back(parse_signature(next(index, name)));
        } else if (name == "--core-max-height") {
            options.core_maximum_height =
                static_cast<SpectralInteger>(std::stoll(next(index, name)));
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
                "unknown projective-core-tail-scan option: " + name);
        }
    }
    if (options.state_paths.empty() ||
        (options.core.empty() && options.core_maximum_height == 0) ||
        (!options.core.empty() && options.core_maximum_height != 0) ||
        options.certificate_path.empty() || options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "projective-core-tail-scan requires states, exactly one core definition, certificate, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveCoreTailScanCli::print_help(std::ostream& out) {
    out << "Projective fixed-core/growing-tail cutoff scan options:\n"
        << "  --state PATH          add one replayable Fourier state; repeatable\n"
        << "  --core-signature A,B,C  add one fixed primitive core ray; repeatable\n"
        << "  --core-max-height H   use every primitive feasible ray with max(a,b,c)<=H\n"
        << "  --certificate PATH    write English JSON scan\n"
        << "  --threads N           parallel projective interaction workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveCoreTailScanCli::run(
    const LocalSldProjectiveCoreTailScanCliOptions& options,
    std::ostream& out) {
    std::vector<SpectralState> states;
    states.reserve(options.state_paths.size());
    for (const std::string& path : options.state_paths) {
        states.push_back(SpectralStateReader::read_tsv(path));
    }
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const std::vector<LocalSldProjectiveCoreSignature> core =
        options.core_maximum_height > 0
        ? ProjectiveCoreFamily::through_maximum_height(
              options.core_maximum_height)
        : options.core;
    const LocalSldProjectiveCoreTailScanReport report =
        LocalSldProjectiveCoreTailScan::analyze(
            dynamics, states, core, options.threads,
            options.exclude_signature_123,
            options.exclude_triple_family);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective core-tail scan certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective core-tail scan rows=" << report.rows.size()
        << " signed_slope="
        << static_cast<double>(report.fitted_open_signed_cutoff_slope)
        << " absolute_slope="
        << static_cast<double>(report.fitted_open_absolute_cutoff_slope)
        << " full_core_signed_slope="
        << static_cast<double>(
               report.fitted_full_active_core_open_signed_cutoff_slope)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_decomposition_exact ? 0 : 2;
}

}  // namespace lemma
