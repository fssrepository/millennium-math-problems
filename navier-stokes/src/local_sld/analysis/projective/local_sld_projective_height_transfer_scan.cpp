#include "local_sld_projective_height_transfer_scan.hpp"

#include "local_sld_projective_height_matrix.hpp"
#include "local_sld_projective_height_tail_summary.hpp"
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
#include <string>

namespace lemma {
namespace {

bool is_power_of_two(SpectralInteger value) {
    return value > 0 && (value & (value - 1)) == 0;
}

SpectralReal fit_slope(
    const std::vector<LocalSldProjectiveHeightTransferRow>& rows,
    bool next_shell) {
    SpectralReal sx = 0.0L;
    SpectralReal sy = 0.0L;
    SpectralReal sxx = 0.0L;
    SpectralReal sxy = 0.0L;
    SpectralReal count = 0.0L;
    for (const auto& row : rows) {
        const SpectralReal value = std::abs(next_shell
            ? row.next_shell_diagonal_power_one
            : row.open_power_one);
        if (value <= 1e-30L) {
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

void write_json(
    const LocalSldProjectiveHeightTransferReport& report,
    const LocalSldProjectiveHeightTransferCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-transfer-v1\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"excludes_signature_123\": "
        << (options.exclude_signature_123 ? "true" : "false")
        << ",\n"
        << "  \"excludes_triple_family\": "
        << (options.exclude_triple_family ? "true" : "false")
        << ",\n"
        << "  \"fitted_open_height_slope\": "
        << static_cast<double>(report.fitted_open_height_slope)
        << ",\n"
        << "  \"fitted_next_shell_diagonal_height_slope\": "
        << static_cast<double>(
               report.fitted_next_shell_diagonal_height_slope)
        << ",\n"
        << "  \"minimum_next_shell_fraction\": "
        << static_cast<double>(report.minimum_next_shell_fraction)
        << ",\n"
        << "  \"maximum_next_shell_fraction\": "
        << static_cast<double>(report.maximum_next_shell_fraction)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"core_maximum_height\": "
            << row.core_maximum_height
            << ", \"cutoff\": " << row.cutoff
            << ", \"open_power_one\": "
            << static_cast<double>(row.open_power_one)
            << ", \"core_tail_power_one\": "
            << static_cast<double>(row.core_tail_power_one)
            << ", \"tail_internal_power_one\": "
            << static_cast<double>(row.tail_internal_power_one)
            << ", \"next_shell_diagonal_power_one\": "
            << static_cast<double>(row.next_shell_diagonal_power_one)
            << ", \"next_shell_fraction_of_open_absolute_sum\": "
            << static_cast<double>(
                   row.next_shell_fraction_of_open_absolute_sum)
            << ", \"open_without_next_shell_diagonal\": "
            << static_cast<double>(
                   row.open_without_next_shell_diagonal)
            << ", \"open_effective_height_pairs\": "
            << static_cast<double>(row.open_effective_height_pairs)
            << ", \"open_signed_alignment\": "
            << static_cast<double>(row.open_signed_alignment)
            << ", \"reconstruction_error\": "
            << static_cast<double>(row.reconstruction_error)
            << ", \"state_path\": \"" << row.state_path << "\"}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_matrix_exact\": "
        << (report.every_matrix_exact ? "true" : "false")
        << ",\n"
        << "  \"finite_scale_transfer_is_not_a_proof\": true,\n"
        << "  \"uniform_weighted_height_bound_proved\": false,\n"
        << "  \"interpretation\": \"the optimized open block follows the first dyadic shell above the fixed core; these finite lower-bound states identify the weighted matrix estimate still required and do not establish a uniform upper bound\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveHeightTransferCliOptions
LocalSldProjectiveHeightTransferCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightTransferCliOptions options;
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
                "unknown projective-height-transfer option: " + name);
        }
    }
    if (options.core_maximum_heights.empty() ||
        options.core_maximum_heights.size() != options.state_paths.size() ||
        options.certificate_path.empty() || options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "projective-height-transfer requires paired heights/states, certificate, and threads 1..256");
    }
    for (const SpectralInteger height : options.core_maximum_heights) {
        if (!is_power_of_two(height)) {
            throw std::invalid_argument(
                "projective-height-transfer heights must be powers of two");
        }
    }
    return options;
}

void LocalSldProjectiveHeightTransferCli::print_help(
    std::ostream& out) {
    out << "Optimized dyadic projective-height transfer options:\n"
        << "  --height H            add a power-of-two core height; repeatable and paired by order\n"
        << "  --state PATH          add its optimized state; repeatable and paired by order\n"
        << "  --certificate PATH    write English JSON transfer scan\n"
        << "  --threads N           parallel direct workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveHeightTransferCli::run(
    const LocalSldProjectiveHeightTransferCliOptions& options,
    std::ostream& out) {
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    LocalSldProjectiveHeightTransferReport report;
    report.every_matrix_exact = true;
    report.minimum_next_shell_fraction =
        std::numeric_limits<SpectralReal>::infinity();
    for (std::size_t index = 0;
         index < options.state_paths.size(); ++index) {
        const SpectralInteger height =
            options.core_maximum_heights[index];
        const SpectralState state = SpectralStateReader::read_tsv(
            options.state_paths[index]);
        const LocalSldProjectiveHeightMatrixReport matrix =
            LocalSldProjectiveHeightMatrix::analyze(
                dynamics, state, options.threads,
                options.exclude_signature_123,
                options.exclude_triple_family);
        const LocalSldProjectiveHeightTailReport tail =
            LocalSldProjectiveHeightTailSummary::summarize(matrix);
        const auto cut = std::find_if(
            tail.rows.begin(), tail.rows.end(),
            [height](const auto& row) {
                return row.core_maximum_height == height;
            });
        if (cut == tail.rows.end()) {
            throw std::invalid_argument(
                "requested height exceeds the active dyadic matrix");
        }
        const int next_shell = cut->last_core_shell + 1;
        const auto diagonal = std::find_if(
            matrix.entries.begin(), matrix.entries.end(),
            [next_shell](const auto& entry) {
                return entry.first_shell == next_shell &&
                    entry.second_shell == next_shell;
            });
        if (diagonal == matrix.entries.end()) {
            throw std::invalid_argument(
                "requested height has no active next dyadic shell");
        }
        LocalSldProjectiveHeightTransferRow row;
        row.core_maximum_height = height;
        row.cutoff = matrix.cutoff;
        row.open_power_one = cut->open_power_one;
        row.core_tail_power_one = cut->core_tail_power_one;
        row.tail_internal_power_one = cut->tail_internal_power_one;
        row.next_shell_diagonal_power_one = diagonal->power_one;
        if (cut->open_absolute_power_one_sum > 0.0L) {
            row.next_shell_fraction_of_open_absolute_sum =
                std::abs(diagonal->power_one) /
                cut->open_absolute_power_one_sum;
        }
        row.open_without_next_shell_diagonal =
            row.open_power_one - row.next_shell_diagonal_power_one;
        row.open_effective_height_pairs =
            cut->open_effective_height_pairs;
        row.open_signed_alignment = cut->open_signed_alignment;
        row.reconstruction_error = std::max(
            matrix.bracket_reconstruction_error,
            cut->reconstruction_error);
        row.state_path = options.state_paths[index];
        report.every_matrix_exact = report.every_matrix_exact &&
            matrix.exact_height_matrix_decomposition &&
            tail.exact_cumulative_decomposition &&
            row.reconstruction_error < 1e-13L;
        report.minimum_next_shell_fraction = std::min(
            report.minimum_next_shell_fraction,
            row.next_shell_fraction_of_open_absolute_sum);
        report.maximum_next_shell_fraction = std::max(
            report.maximum_next_shell_fraction,
            row.next_shell_fraction_of_open_absolute_sum);
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
                "projective height transfer has duplicate heights");
        }
    }
    report.fitted_open_height_slope = fit_slope(report.rows, false);
    report.fitted_next_shell_diagonal_height_slope =
        fit_slope(report.rows, true);

    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective height-transfer certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective height transfer rows=" << report.rows.size()
        << " open_slope="
        << static_cast<double>(report.fitted_open_height_slope)
        << " next_shell_slope="
        << static_cast<double>(
               report.fitted_next_shell_diagonal_height_slope)
        << " next_shell_fraction=["
        << static_cast<double>(report.minimum_next_shell_fraction)
        << ','
        << static_cast<double>(report.maximum_next_shell_fraction)
        << "]\n"
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_matrix_exact ? 0 : 2;
}

}  // namespace lemma
