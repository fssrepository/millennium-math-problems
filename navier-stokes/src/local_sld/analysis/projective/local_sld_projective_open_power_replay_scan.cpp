#include "local_sld_projective_open_power_replay_scan.hpp"

#include "local_sld_projective_open_power_objective.hpp"
#include "projective_core_family.hpp"
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

SpectralReal fit_slope(
    const std::vector<LocalSldProjectiveOpenPowerReplayRow>& rows,
    bool squared) {
    SpectralReal sx = 0.0L;
    SpectralReal sy = 0.0L;
    SpectralReal sxx = 0.0L;
    SpectralReal sxy = 0.0L;
    SpectralReal count = 0.0L;
    for (const auto& row : rows) {
        const SpectralReal value = squared
            ? row.squared_open_power_one
            : row.absolute_open_power_one;
        if (row.core_maximum_height < 1 || value <= 1e-30L) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.core_maximum_height));
        const SpectralReal y = std::log(value);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        count += 1.0L;
    }
    const SpectralReal denominator = count * sxx - sx * sx;
    return count >= 2.0L && std::abs(denominator) > 1e-30L
        ? (count * sxy - sx * sy) / denominator
        : 0.0L;
}

SpectralReal relative_error(
    SpectralReal computed,
    SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

void write_json(
    const LocalSldProjectiveOpenPowerReplayReport& report,
    const LocalSldProjectiveOpenPowerReplayCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-open-power-replay-v1\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"fitted_absolute_height_slope\": "
        << static_cast<double>(report.fitted_absolute_height_slope)
        << ",\n"
        << "  \"fitted_squared_height_slope\": "
        << static_cast<double>(report.fitted_squared_height_slope)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"core_maximum_height\": "
            << row.core_maximum_height
            << ", \"cutoff\": " << row.cutoff
            << ", \"fixed_core_shape_count\": "
            << row.fixed_core_shape_count
            << ", \"absolute_open_power_one\": "
            << static_cast<double>(row.absolute_open_power_one)
            << ", \"squared_open_power_one\": "
            << static_cast<double>(row.squared_open_power_one)
            << ", \"core_tail_power_one\": "
            << static_cast<double>(row.core_tail_power_one)
            << ", \"tail_internal_power_one\": "
            << static_cast<double>(row.tail_internal_power_one)
            << ", \"reconstruction_error\": "
            << static_cast<double>(row.reconstruction_error)
            << ", \"state_path\": \"" << row.state_path << "\"}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_reconstruction_exact\": "
        << (report.every_reconstruction_exact ? "true" : "false")
        << ",\n"
        << "  \"finite_optimized_replay_is_not_a_proof\": true,\n"
        << "  \"uniform_height_tail_bound_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"interpretation\": \"each row uses a state optimized for its own fixed core height; the fitted decay is a finite lower-bound branch, not a uniform upper bound\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveOpenPowerReplayCliOptions
LocalSldProjectiveOpenPowerReplayCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveOpenPowerReplayCliOptions options;
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
                "unknown projective-open-power-replay option: " + name);
        }
    }
    if (options.core_maximum_heights.empty() ||
        options.core_maximum_heights.size() != options.state_paths.size() ||
        options.certificate_path.empty() || options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "projective-open-power-replay requires paired heights/states, certificate, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveOpenPowerReplayCli::print_help(
    std::ostream& out) {
    out << "Optimized projective open-power height replay options:\n"
        << "  --height H            add a fixed core height; repeatable and paired by order\n"
        << "  --state PATH          add its optimized state; repeatable and paired by order\n"
        << "  --certificate PATH    write English JSON replay scan\n"
        << "  --threads N           parallel direct workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveOpenPowerReplayCli::run(
    const LocalSldProjectiveOpenPowerReplayCliOptions& options,
    std::ostream& out) {
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const TriadSelection selection = selection_for(
        options.exclude_signature_123,
        options.exclude_triple_family);
    LocalSldProjectiveOpenPowerReplayReport report;
    report.every_reconstruction_exact = true;
    for (std::size_t index = 0;
         index < options.state_paths.size(); ++index) {
        const SpectralInteger height =
            options.core_maximum_heights[index];
        const SpectralState state = SpectralStateReader::read_tsv(
            options.state_paths[index]);
        const LocalSldProjectiveOpenPowerObjectiveValue value =
            LocalSldProjectiveOpenPowerObjective(
                dynamics, selection, height, options.threads)
                .evaluate(state);
        const auto core = ProjectiveCoreFamily::through_maximum_height(
            height);
        const LocalSldProjectiveCoreTailReport decomposition =
            LocalSldProjectiveCoreTailLedger::analyze(
                dynamics, state, core, options.threads,
                options.exclude_signature_123,
                options.exclude_triple_family);
        LocalSldProjectiveOpenPowerReplayRow row;
        row.core_maximum_height = height;
        row.cutoff = SpectralStateOps::cutoff(state);
        row.fixed_core_shape_count = core.size();
        row.absolute_open_power_one = value.absolute_open_power_one;
        row.squared_open_power_one = value.squared_open_power_one;
        row.core_tail_power_one = decomposition.core_tail.power_one;
        row.tail_internal_power_one = decomposition.tail.power_one;
        row.reconstruction_error = relative_error(
            std::abs(row.core_tail_power_one +
                     row.tail_internal_power_one),
            row.absolute_open_power_one);
        row.state_path = options.state_paths[index];
        report.every_reconstruction_exact =
            report.every_reconstruction_exact &&
            decomposition.exact_core_tail_decomposition &&
            row.reconstruction_error < 1e-13L;
        report.rows.push_back(row);
    }
    std::sort(report.rows.begin(), report.rows.end(),
        [](const auto& left, const auto& right) {
            return left.core_maximum_height <
                right.core_maximum_height;
        });
    for (std::size_t index = 1; index < report.rows.size(); ++index) {
        if (report.rows[index - 1].core_maximum_height ==
            report.rows[index].core_maximum_height) {
            throw std::invalid_argument(
                "projective open-power replay has duplicate heights");
        }
    }
    report.fitted_absolute_height_slope = fit_slope(
        report.rows, false);
    report.fitted_squared_height_slope = fit_slope(
        report.rows, true);

    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective open-power replay certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective open-power replay rows=" << report.rows.size()
        << " absolute_height_slope="
        << static_cast<double>(report.fitted_absolute_height_slope)
        << " squared_height_slope="
        << static_cast<double>(report.fitted_squared_height_slope)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_reconstruction_exact ? 0 : 2;
}

}  // namespace lemma
