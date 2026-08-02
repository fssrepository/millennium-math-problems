#include "local_sld_projective_normalization_cauchy_scan.hpp"

#include "local_sld_projective_normalization_cauchy_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal fitted_log_slope(
    const std::vector<LocalSldProjectiveNormalizationCauchyScanRow>& rows) {
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    std::size_t count = 0;
    for (const auto& row : rows) {
        if (row.core_maximum_height <= 0 ||
            !(row.cauchy_bound_power_one > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.core_maximum_height));
        const SpectralReal y = std::log(row.cauchy_bound_power_one);
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
    const SpectralReal divisor = n * sum_xx - sum_x * sum_x;
    return std::abs(divisor) > 1e-30L
        ? (n * sum_xy - sum_x * sum_y) / divisor
        : 0.0L;
}

void write_json(
    const LocalSldProjectiveNormalizationCauchyScanReport& report,
    const LocalSldProjectiveNormalizationCauchyScanOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-cauchy-height-scan-v1\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"fitted_height_slope\": "
        << static_cast<double>(report.fitted_height_slope) << ",\n"
        << "  \"maximum_bound\": "
        << static_cast<double>(report.maximum_bound) << ",\n"
        << "  \"minimum_bound\": "
        << static_cast<double>(report.minimum_bound) << ",\n"
        << "  \"fitted_quarter_compensated_slope\": "
        << static_cast<double>(
               report.fitted_quarter_compensated_slope) << ",\n"
        << "  \"maximum_quarter_compensated_bound\": "
        << static_cast<double>(
               report.maximum_quarter_compensated_bound) << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"core_maximum_height\": "
            << row.core_maximum_height
            << ", \"cutoff\": " << row.cutoff
            << ", \"full_stretching\": "
            << static_cast<double>(row.full_stretching)
            << ", \"enstrophy\": "
            << static_cast<double>(row.enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(row.palinstrophy)
            << ", \"selected_aggregate_h1_norm2\": "
            << static_cast<double>(row.selected_aggregate_h1_norm2)
            << ", \"tail_aggregate_h2_norm2\": "
            << static_cast<double>(row.tail_aggregate_h2_norm2)
            << ", \"cauchy_bound_power_one\": "
            << static_cast<double>(row.cauchy_bound_power_one)
            << ", \"squared_cauchy_bound_power_one\": "
            << static_cast<double>(row.squared_cauchy_bound_power_one)
            << ", \"quarter_compensated_bound\": "
            << static_cast<double>(row.quarter_compensated_bound)
            << ", \"half_compensated_squared_bound\": "
            << static_cast<double>(row.half_compensated_squared_bound)
            << ", \"numerator\": "
            << static_cast<double>(row.numerator)
            << ", \"denominator\": "
            << static_cast<double>(row.denominator)
            << ", \"state_path\": \"" << row.state_path
            << "\", \"finite\": " << (row.finite ? "true" : "false")
            << '}' << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_row_finite\": "
        << (report.every_row_finite ? "true" : "false") << ",\n"
        << "  \"finite_scan_is_not_a_proof\": true,\n"
        << "  \"uniform_height_decay_proved\": false,\n"
        << "  \"interpretation\": \"finite optimized height scaling of the exact canonical Cauchy majorant; a negative fitted slope is a lemma candidate, not a uniform estimate\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveNormalizationCauchyScanOptions
LocalSldProjectiveNormalizationCauchyScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationCauchyScanOptions options;
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
                "unknown normalization-Cauchy-scan option: " + name);
        }
    }
    if (options.core_maximum_heights.size() < 2 ||
        options.core_maximum_heights.size() != options.state_paths.size() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.threads < 1 || options.threads > 256 ||
        std::any_of(
            options.core_maximum_heights.begin(),
            options.core_maximum_heights.end(),
            [](SpectralInteger height) { return height < 1; })) {
        throw std::invalid_argument(
            "normalization-Cauchy-scan requires at least two paired positive heights/states, a certificate, valid selection, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveNormalizationCauchyScanCli::print_help(
    std::ostream& out) {
    out << "Projective normalization Cauchy-majorant height scan options:\n"
        << "  --height H          add a core height; repeatable and paired by order\n"
        << "  --state PATH        add its optimized state; repeatable and paired by order\n"
        << "  --selection NAME    local SLD triad selection\n"
        << "  --threads N         parallel direct workers\n"
        << "  --certificate PATH  write English JSON scan\n";
}

int LocalSldProjectiveNormalizationCauchyScanCli::run(
    const LocalSldProjectiveNormalizationCauchyScanOptions& options,
    std::ostream& out) {
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    LocalSldProjectiveNormalizationCauchyScanReport report;
    report.every_row_finite = true;
    report.minimum_bound = std::numeric_limits<SpectralReal>::infinity();
    for (std::size_t index = 0; index < options.state_paths.size(); ++index) {
        const SpectralState state = SpectralStateReader::read_tsv(
            options.state_paths[index]);
        const auto value = LocalSldProjectiveNormalizationCauchyObjective(
            dynamics,
            LocalSldTriadSelection::parse(options.selection),
            options.core_maximum_heights[index],
            options.threads).evaluate(state);
        LocalSldProjectiveNormalizationCauchyScanRow row;
        row.core_maximum_height = options.core_maximum_heights[index];
        row.cutoff = SpectralStateOps::cutoff(state);
        row.full_stretching = value.full_stretching;
        row.enstrophy = value.enstrophy;
        row.palinstrophy = value.palinstrophy;
        row.selected_aggregate_h1_norm2 =
            value.selected_aggregate_h1_norm2;
        row.tail_aggregate_h2_norm2 = value.tail_aggregate_h2_norm2;
        row.cauchy_bound_power_one = value.cauchy_bound_power_one;
        row.squared_cauchy_bound_power_one =
            value.squared_cauchy_bound_power_one;
        row.quarter_compensated_bound = std::pow(
            static_cast<SpectralReal>(row.core_maximum_height), 0.25L) *
            row.cauchy_bound_power_one;
        row.half_compensated_squared_bound = std::sqrt(
            static_cast<SpectralReal>(row.core_maximum_height)) *
            row.squared_cauchy_bound_power_one;
        row.numerator = 1.5L * std::abs(row.full_stretching) * std::sqrt(
            std::max(
                row.enstrophy * row.selected_aggregate_h1_norm2 *
                    row.palinstrophy * row.tail_aggregate_h2_norm2,
                0.0L));
        row.denominator =
            row.enstrophy * row.enstrophy *
            row.palinstrophy * row.palinstrophy * row.palinstrophy;
        row.state_path = options.state_paths[index];
        row.finite = value.finite && std::isfinite(row.numerator) &&
            std::isfinite(row.denominator);
        report.every_row_finite = report.every_row_finite && row.finite;
        report.maximum_bound = std::max(
            report.maximum_bound, row.cauchy_bound_power_one);
        report.minimum_bound = std::min(
            report.minimum_bound, row.cauchy_bound_power_one);
        report.maximum_quarter_compensated_bound = std::max(
            report.maximum_quarter_compensated_bound,
            row.quarter_compensated_bound);
        report.rows.push_back(std::move(row));
    }
    report.fitted_height_slope = fitted_log_slope(report.rows);
    report.fitted_quarter_compensated_slope =
        report.fitted_height_slope + 0.25L;
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write normalization-Cauchy scan certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective normalization Cauchy scan rows="
        << report.rows.size()
        << " height_slope="
        << static_cast<double>(report.fitted_height_slope)
        << " compensated_slope="
        << static_cast<double>(
               report.fitted_quarter_compensated_slope)
        << " max_H1/4_bound="
        << static_cast<double>(
               report.maximum_quarter_compensated_bound)
        << " bound=[" << static_cast<double>(report.minimum_bound)
        << ',' << static_cast<double>(report.maximum_bound) << "]\n"
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_row_finite ? 0 : 2;
}

}  // namespace lemma
