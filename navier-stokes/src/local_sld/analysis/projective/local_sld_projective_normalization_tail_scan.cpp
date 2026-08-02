#include "local_sld_projective_normalization_tail_scan.hpp"

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

enum class Series {
    open,
    core_stretching_tail_cross,
    tail_stretching_core_cross,
    tail_stretching_tail_cross,
    joint_cross_tail_cauchy_bound,
};

SpectralReal series_value(
    const LocalSldProjectiveNormalizationTailRow& row,
    Series series) {
    switch (series) {
        case Series::open:
            return row.open_palinstrophy_normalization_power_one;
        case Series::core_stretching_tail_cross:
            return row.core_stretching_tail_cross_power_one;
        case Series::tail_stretching_core_cross:
            return row.tail_stretching_core_cross_power_one;
        case Series::tail_stretching_tail_cross:
            return row.tail_stretching_tail_cross_power_one;
        case Series::joint_cross_tail_cauchy_bound:
            return row.joint_cross_tail_cauchy_bound;
    }
    return 0.0L;
}

SpectralReal fit_slope(
    const std::vector<LocalSldProjectiveNormalizationTailRow>& rows,
    Series series) {
    SpectralReal sx = 0.0L;
    SpectralReal sy = 0.0L;
    SpectralReal sxx = 0.0L;
    SpectralReal sxy = 0.0L;
    SpectralReal count = 0.0L;
    for (const auto& row : rows) {
        const SpectralReal value = std::abs(series_value(row, series));
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

std::string dominant_channel(
    const LocalSldProjectiveNormalizationTailRow& row) {
    const SpectralReal first = std::abs(
        row.core_stretching_tail_cross_power_one);
    const SpectralReal second = std::abs(
        row.tail_stretching_core_cross_power_one);
    const SpectralReal third = std::abs(
        row.tail_stretching_tail_cross_power_one);
    if (first >= second && first >= third) {
        return "core-stretching-times-tail-cross";
    }
    if (second >= third) {
        return "tail-stretching-times-core-cross";
    }
    return "tail-stretching-times-tail-cross";
}

void write_json(
    const LocalSldProjectiveNormalizationTailReport& report,
    const LocalSldProjectiveNormalizationTailCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-tail-scan-v1\",\n"
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
        << "  \"fitted_core_stretching_tail_cross_height_slope\": "
        << static_cast<double>(
               report.fitted_core_stretching_tail_cross_height_slope)
        << ",\n"
        << "  \"fitted_tail_stretching_core_cross_height_slope\": "
        << static_cast<double>(
               report.fitted_tail_stretching_core_cross_height_slope)
        << ",\n"
        << "  \"fitted_tail_stretching_tail_cross_height_slope\": "
        << static_cast<double>(
               report.fitted_tail_stretching_tail_cross_height_slope)
        << ",\n"
        << "  \"fitted_joint_cross_tail_cauchy_bound_height_slope\": "
        << static_cast<double>(
               report.fitted_joint_cross_tail_cauchy_bound_height_slope)
        << ",\n"
        << "  \"maximum_factorization_error\": "
        << static_cast<double>(report.maximum_factorization_error)
        << ",\n"
        << "  \"minimum_joint_cross_tail_alignment\": "
        << static_cast<double>(report.minimum_joint_cross_tail_alignment)
        << ",\n"
        << "  \"maximum_joint_cross_tail_alignment\": "
        << static_cast<double>(report.maximum_joint_cross_tail_alignment)
        << ",\n"
        << "  \"minimum_joint_cross_tail_cauchy_ratio\": "
        << static_cast<double>(report.minimum_joint_cross_tail_cauchy_ratio)
        << ",\n"
        << "  \"maximum_joint_cross_tail_cauchy_ratio\": "
        << static_cast<double>(report.maximum_joint_cross_tail_cauchy_ratio)
        << ",\n"
        << "  \"maximum_individual_cauchy_ratio\": "
        << static_cast<double>(report.maximum_individual_cauchy_ratio)
        << ",\n"
        << "  \"dominance_counts\": {\n"
        << "    \"core_stretching_times_tail_cross\": "
        << report.core_stretching_tail_cross_dominance_count << ",\n"
        << "    \"tail_stretching_times_core_cross\": "
        << report.tail_stretching_core_cross_dominance_count << ",\n"
        << "    \"tail_stretching_times_tail_cross\": "
        << report.tail_stretching_tail_cross_dominance_count << "\n"
        << "  },\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"core_maximum_height\": "
            << row.core_maximum_height
            << ", \"cutoff\": " << row.cutoff
            << ", \"open_palinstrophy_normalization_power_one\": "
            << static_cast<double>(
                   row.open_palinstrophy_normalization_power_one)
            << ", \"core_stretching_tail_cross_power_one\": "
            << static_cast<double>(
                   row.core_stretching_tail_cross_power_one)
            << ", \"tail_stretching_core_cross_power_one\": "
            << static_cast<double>(
                   row.tail_stretching_core_cross_power_one)
            << ", \"tail_stretching_tail_cross_power_one\": "
            << static_cast<double>(
                   row.tail_stretching_tail_cross_power_one)
            << ", \"joint_cross_tail_absolute_envelope\": "
            << static_cast<double>(
                   row.joint_cross_tail_absolute_envelope)
            << ", \"joint_cross_tail_alignment\": "
            << static_cast<double>(row.joint_cross_tail_alignment)
            << ", \"core_stretching_tail_cross_fraction\": "
            << static_cast<double>(
                   row.core_stretching_tail_cross_fraction)
            << ", \"tail_stretching_core_cross_fraction\": "
            << static_cast<double>(
                   row.tail_stretching_core_cross_fraction)
            << ", \"tail_stretching_tail_cross_fraction\": "
            << static_cast<double>(
                   row.tail_stretching_tail_cross_fraction)
            << ", \"core_stretching_tail_cross_cauchy_bound\": "
            << static_cast<double>(
                   row.core_stretching_tail_cross_cauchy_bound)
            << ", \"tail_stretching_core_cross_cauchy_bound\": "
            << static_cast<double>(
                   row.tail_stretching_core_cross_cauchy_bound)
            << ", \"tail_stretching_tail_cross_cauchy_bound\": "
            << static_cast<double>(
                   row.tail_stretching_tail_cross_cauchy_bound)
            << ", \"core_stretching_tail_cross_cauchy_ratio\": "
            << static_cast<double>(
                   row.core_stretching_tail_cross_cauchy_ratio)
            << ", \"tail_stretching_core_cross_cauchy_ratio\": "
            << static_cast<double>(
                   row.tail_stretching_core_cross_cauchy_ratio)
            << ", \"tail_stretching_tail_cross_cauchy_ratio\": "
            << static_cast<double>(
                   row.tail_stretching_tail_cross_cauchy_ratio)
            << ", \"joint_cross_tail_cauchy_bound\": "
            << static_cast<double>(row.joint_cross_tail_cauchy_bound)
            << ", \"joint_cross_tail_cauchy_ratio\": "
            << static_cast<double>(row.joint_cross_tail_cauchy_ratio)
            << ", \"core_stretching_h1_alignment\": "
            << static_cast<double>(row.core_stretching_h1_alignment)
            << ", \"tail_stretching_h1_alignment\": "
            << static_cast<double>(row.tail_stretching_h1_alignment)
            << ", \"core_palinstrophy_cross_h2_alignment\": "
            << static_cast<double>(
                   row.core_palinstrophy_cross_h2_alignment)
            << ", \"tail_palinstrophy_cross_h2_alignment\": "
            << static_cast<double>(
                   row.tail_palinstrophy_cross_h2_alignment)
            << ", \"core_stretching\": "
            << static_cast<double>(row.core_stretching)
            << ", \"tail_stretching\": "
            << static_cast<double>(row.tail_stretching)
            << ", \"core_palinstrophy_cross\": "
            << static_cast<double>(row.core_palinstrophy_cross)
            << ", \"tail_palinstrophy_cross\": "
            << static_cast<double>(row.tail_palinstrophy_cross)
            << ", \"factorization_error\": "
            << static_cast<double>(row.factorization_error)
            << ", \"dominant_channel\": \""
            << row.dominant_channel
            << "\", \"state_path\": \"" << row.state_path << "\"}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"dominant_channel_switch_observed\": "
        << (report.dominant_channel_switch_observed ? "true" : "false")
        << ",\n"
        << "  \"every_matrix_exact\": "
        << (report.every_matrix_exact ? "true" : "false")
        << ",\n"
        << "  \"finite_optimized_scan_is_not_a_proof\": true,\n"
        << "  \"uniform_joint_cross_tail_bound_proved\": false,\n"
        << "  \"candidate_lemma\": \"uniformly control the joint normalized absolute sum of s_core*t_tail, s_tail*t_core, and s_tail*t_tail; no single channel can be discarded on the optimized finite family\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveNormalizationTailCliOptions
LocalSldProjectiveNormalizationTailScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationTailCliOptions options;
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
                "unknown normalization-tail-scan option: " + name);
        }
    }
    if (options.core_maximum_heights.empty() ||
        options.core_maximum_heights.size() != options.state_paths.size() ||
        options.certificate_path.empty() || options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "normalization-tail-scan requires paired heights/states, certificate, and threads 1..256");
    }
    for (const SpectralInteger height : options.core_maximum_heights) {
        if (!is_power_of_two(height)) {
            throw std::invalid_argument(
                "normalization-tail-scan heights must be powers of two");
        }
    }
    return options;
}

void LocalSldProjectiveNormalizationTailScanCli::print_help(
    std::ostream& out) {
    out << "Projective normalization core/tail factor scan options:\n"
        << "  --height H            add a power-of-two core height; repeatable and paired by order\n"
        << "  --state PATH          add its optimized state; repeatable and paired by order\n"
        << "  --certificate PATH    write English JSON factor scan\n"
        << "  --threads N           parallel direct workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveNormalizationTailScanCli::run(
    const LocalSldProjectiveNormalizationTailCliOptions& options,
    std::ostream& out) {
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    LocalSldProjectiveNormalizationTailReport report;
    report.every_matrix_exact = true;
    report.minimum_joint_cross_tail_alignment =
        std::numeric_limits<SpectralReal>::infinity();
    report.minimum_joint_cross_tail_cauchy_ratio =
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
        LocalSldProjectiveNormalizationTailRow row;
        row.core_maximum_height = height;
        row.cutoff = matrix.cutoff;
        row.open_palinstrophy_normalization_power_one =
            cut->open_palinstrophy_normalization_power_one;
        row.core_stretching_tail_cross_power_one =
            cut->core_stretching_tail_cross_power_one;
        row.tail_stretching_core_cross_power_one =
            cut->tail_stretching_core_cross_power_one;
        row.tail_stretching_tail_cross_power_one =
            cut->tail_stretching_tail_cross_power_one;
        row.joint_cross_tail_absolute_envelope =
            std::abs(row.core_stretching_tail_cross_power_one) +
            std::abs(row.tail_stretching_core_cross_power_one) +
            std::abs(row.tail_stretching_tail_cross_power_one);
        if (row.joint_cross_tail_absolute_envelope > 0.0L) {
            const SpectralReal inverse =
                1.0L / row.joint_cross_tail_absolute_envelope;
            row.joint_cross_tail_alignment =
                std::abs(
                    row.open_palinstrophy_normalization_power_one) *
                inverse;
            row.core_stretching_tail_cross_fraction =
                std::abs(
                    row.core_stretching_tail_cross_power_one) * inverse;
            row.tail_stretching_core_cross_fraction =
                std::abs(
                    row.tail_stretching_core_cross_power_one) * inverse;
        row.tail_stretching_tail_cross_fraction =
                std::abs(
                    row.tail_stretching_tail_cross_power_one) * inverse;
        }
        row.core_stretching_tail_cross_cauchy_bound =
            cut->core_stretching_tail_cross_cauchy_bound;
        row.tail_stretching_core_cross_cauchy_bound =
            cut->tail_stretching_core_cross_cauchy_bound;
        row.tail_stretching_tail_cross_cauchy_bound =
            cut->tail_stretching_tail_cross_cauchy_bound;
        row.core_stretching_tail_cross_cauchy_ratio =
            cut->core_stretching_tail_cross_cauchy_ratio;
        row.tail_stretching_core_cross_cauchy_ratio =
            cut->tail_stretching_core_cross_cauchy_ratio;
        row.tail_stretching_tail_cross_cauchy_ratio =
            cut->tail_stretching_tail_cross_cauchy_ratio;
        row.joint_cross_tail_cauchy_bound =
            cut->joint_cross_tail_cauchy_bound;
        row.joint_cross_tail_cauchy_ratio =
            cut->joint_cross_tail_cauchy_ratio;
        row.core_stretching_h1_alignment =
            cut->core_stretching_h1_alignment;
        row.tail_stretching_h1_alignment =
            cut->tail_stretching_h1_alignment;
        row.core_palinstrophy_cross_h2_alignment =
            cut->core_palinstrophy_cross_h2_alignment;
        row.tail_palinstrophy_cross_h2_alignment =
            cut->tail_palinstrophy_cross_h2_alignment;
        row.core_stretching = cut->core_stretching;
        row.tail_stretching = cut->tail_stretching;
        row.core_palinstrophy_cross = cut->core_palinstrophy_cross;
        row.tail_palinstrophy_cross = cut->tail_palinstrophy_cross;
        row.factorization_error = cut->palinstrophy_factorization_error;
        row.dominant_channel = dominant_channel(row);
        row.state_path = options.state_paths[index];
        report.maximum_factorization_error = std::max(
            report.maximum_factorization_error,
            row.factorization_error);
        report.minimum_joint_cross_tail_alignment = std::min(
            report.minimum_joint_cross_tail_alignment,
            row.joint_cross_tail_alignment);
        report.maximum_joint_cross_tail_alignment = std::max(
            report.maximum_joint_cross_tail_alignment,
            row.joint_cross_tail_alignment);
        report.minimum_joint_cross_tail_cauchy_ratio = std::min(
            report.minimum_joint_cross_tail_cauchy_ratio,
            row.joint_cross_tail_cauchy_ratio);
        report.maximum_joint_cross_tail_cauchy_ratio = std::max(
            report.maximum_joint_cross_tail_cauchy_ratio,
            row.joint_cross_tail_cauchy_ratio);
        report.maximum_individual_cauchy_ratio = std::max({
            report.maximum_individual_cauchy_ratio,
            row.core_stretching_tail_cross_cauchy_ratio,
            row.tail_stretching_core_cross_cauchy_ratio,
            row.tail_stretching_tail_cross_cauchy_ratio});
        if (row.dominant_channel ==
            "core-stretching-times-tail-cross") {
            ++report.core_stretching_tail_cross_dominance_count;
        } else if (row.dominant_channel ==
                   "tail-stretching-times-core-cross") {
            ++report.tail_stretching_core_cross_dominance_count;
        } else {
            ++report.tail_stretching_tail_cross_dominance_count;
        }
        report.every_matrix_exact = report.every_matrix_exact &&
            matrix.exact_height_matrix_decomposition &&
            tail.exact_cumulative_decomposition &&
            row.factorization_error < 1e-13L;
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
                "normalization tail scan has duplicate heights");
        }
        report.dominant_channel_switch_observed =
            report.dominant_channel_switch_observed ||
            report.rows[index - 1].dominant_channel !=
                report.rows[index].dominant_channel;
    }
    report.fitted_open_height_slope = fit_slope(
        report.rows, Series::open);
    report.fitted_core_stretching_tail_cross_height_slope = fit_slope(
        report.rows, Series::core_stretching_tail_cross);
    report.fitted_tail_stretching_core_cross_height_slope = fit_slope(
        report.rows, Series::tail_stretching_core_cross);
    report.fitted_tail_stretching_tail_cross_height_slope = fit_slope(
        report.rows, Series::tail_stretching_tail_cross);
    report.fitted_joint_cross_tail_cauchy_bound_height_slope = fit_slope(
        report.rows, Series::joint_cross_tail_cauchy_bound);

    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write normalization tail-scan certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective normalization tail scan rows="
        << report.rows.size()
        << " open_slope="
        << static_cast<double>(report.fitted_open_height_slope)
        << " factor_slopes=["
        << static_cast<double>(
               report.fitted_core_stretching_tail_cross_height_slope)
        << ','
        << static_cast<double>(
               report.fitted_tail_stretching_core_cross_height_slope)
        << ','
        << static_cast<double>(
               report.fitted_tail_stretching_tail_cross_height_slope)
        << "] cauchy_bound_slope="
        << static_cast<double>(
               report.fitted_joint_cross_tail_cauchy_bound_height_slope)
        << " cauchy_ratio=["
        << static_cast<double>(
               report.minimum_joint_cross_tail_cauchy_ratio)
        << ','
        << static_cast<double>(
               report.maximum_joint_cross_tail_cauchy_ratio)
        << "] switch="
        << (report.dominant_channel_switch_observed ? "yes" : "no")
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_matrix_exact ? 0 : 2;
}

}  // namespace lemma
