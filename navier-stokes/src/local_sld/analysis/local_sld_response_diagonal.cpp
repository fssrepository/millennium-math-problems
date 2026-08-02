#include "local_sld_response_diagonal.hpp"

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

void write_real_array(
    std::ostream& output,
    const std::vector<SpectralReal>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        output << static_cast<double>(values[index]);
        if (index + 1 != values.size()) {
            output << ", ";
        }
    }
    output << ']';
}

void write_int_array(
    std::ostream& output,
    const std::vector<int>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        output << values[index];
        if (index + 1 != values.size()) {
            output << ", ";
        }
    }
    output << ']';
}

void write_certificate(
    const LocalSldResponseDiagonalReport& report,
    const LocalSldResponseDiagonalOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write response-diagonal certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-response-diagonal-v2\",\n"
        << "  \"definition\": \"at cutoff K retain only response orders 0 through K, before the recursion requires shell K+1\",\n"
        << "  \"weighted_l1_definition\": \"A_r(K)=sum_{n=0}^K r^n |<u,b_n>|/sqrt(E)\",\n"
        << "  \"quadratic_identity\": \"sum r^n sum_{i+j=n-1} a_i a_j <= r A_r(K)^2\",\n"
        << "  \"maximum_depth\": " << report.maximum_depth << ",\n"
        << "  \"radius\": " << static_cast<double>(report.radius)
        << ",\n"
        << "  \"maximum_safe_weighted_l1\": "
        << static_cast<double>(report.maximum_safe_weighted_l1) << ",\n"
        << "  \"maximum_common_coefficient_l2_difference\": "
        << static_cast<double>(
               report.maximum_common_coefficient_l2_difference)
        << ",\n"
        << "  \"maximum_quadratic_majorant_ratio\": "
        << static_cast<double>(report.maximum_quadratic_majorant_ratio)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const LocalSldResponseDiagonalRow& row = report.rows[index];
        output << "    {\"state_path\": \"" << row.state_path
            << "\", \"cutoff\": " << row.cutoff
            << ", \"safe_order_count\": " << row.safe_order_count
            << ", \"first_truncated_order\": "
            << row.first_truncated_order
            << ", \"safe_shells_inside_cutoff\": "
            << (row.safe_shells_inside_cutoff ? "true" : "false")
            << ", \"safe_projection_energy\": "
            << static_cast<double>(row.safe_projection_energy)
            << ", \"safe_weighted_l1\": "
            << static_cast<double>(row.safe_weighted_l1)
            << ", \"safe_weighted_tail_after_cubic\": "
            << static_cast<double>(row.safe_weighted_tail_after_cubic)
            << ", \"last_safe_coefficient\": "
            << static_cast<double>(row.last_safe_coefficient)
            << ", \"first_truncated_coefficient\": "
            << static_cast<double>(row.first_truncated_coefficient)
            << ", \"quadratic_majorant\": "
            << static_cast<double>(row.quadratic_majorant)
            << ", \"quadratic_majorant_bound\": "
            << static_cast<double>(row.quadratic_majorant_bound)
            << ", \"quadratic_majorant_ratio\": "
            << static_cast<double>(row.quadratic_majorant_ratio)
            << ", \"common_coefficient_l2_difference\": "
            << static_cast<double>(
                   row.common_coefficient_l2_difference)
            << ", \"common_weighted_l1_difference\": "
            << static_cast<double>(
                   row.common_weighted_l1_difference)
            << ", \"safe_coefficients\": ";
        write_real_array(output, row.safe_coefficients);
        output << ", \"safe_highest_shells\": ";
        write_int_array(output, row.safe_highest_shells);
        output << ", \"safe_lowest_shells\": ";
        write_int_array(output, row.safe_lowest_shells);
        output << ", \"safe_highest_shell_energy_fractions\": ";
        write_real_array(
            output, row.safe_highest_shell_energy_fractions);
        output << ", \"safe_lower_half_shell_energy_fractions\": ";
        write_real_array(
            output, row.safe_lower_half_shell_energy_fractions);
        output << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_cutoff_diagonal_is_not_a_proof\": true,\n"
        << "  \"remaining_requirement\": \"prove a cutoff-uniform A_r bound and a compatible bilinear response norm bound for some r greater than one\"\n"
        << "}\n";
}

}  // namespace

LocalSldResponseDiagonalReport LocalSldResponseDiagonal::analyze(
    const SpectralDynamics& dynamics,
    const std::vector<std::pair<std::string, SpectralState>>& states,
    int maximum_depth,
    SpectralReal radius) {
    if (states.empty() || maximum_depth < 2 || maximum_depth > 16 ||
        !(radius >= 1.0L) || !std::isfinite(radius)) {
        throw std::invalid_argument(
            "invalid response-diagonal configuration");
    }
    LocalSldResponseDiagonalReport report;
    report.maximum_depth = maximum_depth;
    report.radius = radius;
    for (const auto& [path, state] : states) {
        const int cutoff = SpectralStateOps::cutoff(state);
        if (cutoff < 2) {
            throw std::invalid_argument(
                "response diagonal requires cutoff at least two");
        }
        const int construction_depth = std::min(
            maximum_depth, cutoff + 2);
        const LocalSldResponseHierarchyReport hierarchy =
            LocalSldResponseHierarchy::analyze(
                dynamics, state, construction_depth, false, false);
        LocalSldResponseDiagonalRow row;
        row.state_path = path;
        row.cutoff = cutoff;
        row.safe_order_count = std::min(
            static_cast<int>(hierarchy.rows.size()), cutoff + 1);
        row.first_truncated_order =
            static_cast<int>(hierarchy.rows.size()) > row.safe_order_count
            ? row.safe_order_count
            : -1;
        const SpectralReal inverse_energy_norm = 1.0L / std::sqrt(
            std::max(hierarchy.reference_energy, 1e-30L));
        row.safe_coefficients.reserve(
            static_cast<std::size_t>(row.safe_order_count));
        row.safe_highest_shells.reserve(
            static_cast<std::size_t>(row.safe_order_count));
        row.safe_lowest_shells.reserve(
            static_cast<std::size_t>(row.safe_order_count));
        row.safe_highest_shell_energy_fractions.reserve(
            static_cast<std::size_t>(row.safe_order_count));
        row.safe_lower_half_shell_energy_fractions.reserve(
            static_cast<std::size_t>(row.safe_order_count));
        for (int order = 0; order < row.safe_order_count; ++order) {
            const LocalSldResponseHierarchyRow& entry =
                hierarchy.rows[static_cast<std::size_t>(order)];
            const SpectralReal coefficient =
                std::abs(entry.coefficient) * inverse_energy_norm;
            const SpectralReal weight = std::pow(radius, order);
            row.safe_coefficients.push_back(coefficient);
            row.safe_highest_shells.push_back(
                entry.highest_active_shell);
            row.safe_lowest_shells.push_back(
                entry.lowest_active_shell);
            row.safe_highest_shell_energy_fractions.push_back(
                entry.highest_shell_energy_fraction);
            row.safe_lower_half_shell_energy_fractions.push_back(
                entry.lower_half_shell_energy_fraction);
            row.safe_shells_inside_cutoff =
                row.safe_shells_inside_cutoff &&
                entry.highest_active_shell <= cutoff;
            row.safe_projection_energy += entry.coefficient_energy;
            row.safe_weighted_l1 += weight * coefficient;
            if (order >= 3) {
                row.safe_weighted_tail_after_cubic +=
                    weight * coefficient;
            }
        }
        row.last_safe_coefficient = row.safe_coefficients.back();
        if (row.first_truncated_order >= 0) {
            row.first_truncated_coefficient = std::abs(
                hierarchy.rows[static_cast<std::size_t>(
                    row.first_truncated_order)].coefficient) *
                inverse_energy_norm;
        }
        for (int order = 1; order < row.safe_order_count; ++order) {
            SpectralReal convolution = 0.0L;
            for (int left = 0; left < order; ++left) {
                const int right = order - 1 - left;
                convolution += row.safe_coefficients[
                    static_cast<std::size_t>(left)] *
                    row.safe_coefficients[
                        static_cast<std::size_t>(right)];
            }
            row.quadratic_majorant +=
                std::pow(radius, order) * convolution;
        }
        row.quadratic_majorant_bound =
            radius * row.safe_weighted_l1 * row.safe_weighted_l1;
        row.quadratic_majorant_ratio = row.quadratic_majorant_bound > 0.0L
            ? row.quadratic_majorant / row.quadratic_majorant_bound
            : 0.0L;
        if (!report.rows.empty()) {
            const LocalSldResponseDiagonalRow& previous =
                report.rows.back();
            const std::size_t common = std::min(
                previous.safe_coefficients.size(),
                row.safe_coefficients.size());
            SpectralReal difference2 = 0.0L;
            for (std::size_t order = 0; order < common; ++order) {
                const SpectralReal difference =
                    row.safe_coefficients[order] -
                    previous.safe_coefficients[order];
                difference2 += difference * difference;
                row.common_weighted_l1_difference +=
                    std::pow(radius, static_cast<int>(order)) *
                    std::abs(difference);
            }
            row.common_coefficient_l2_difference =
                std::sqrt(difference2);
        }
        report.maximum_safe_weighted_l1 = std::max(
            report.maximum_safe_weighted_l1,
            row.safe_weighted_l1);
        report.maximum_common_coefficient_l2_difference = std::max(
            report.maximum_common_coefficient_l2_difference,
            row.common_coefficient_l2_difference);
        report.maximum_quadratic_majorant_ratio = std::max(
            report.maximum_quadratic_majorant_ratio,
            row.quadratic_majorant_ratio);
        report.rows.push_back(std::move(row));
    }
    return report;
}

LocalSldResponseDiagonalOptions LocalSldResponseDiagonalCli::parse(
    int argc, char** argv, int first) {
    LocalSldResponseDiagonalOptions options;
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
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--max-depth") {
            options.maximum_depth = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--radius") {
            options.radius = std::stold(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown local-sld-response-diagonal option: " + name);
        }
    }
    if (options.state_paths.empty() || options.certificate_path.empty() ||
        options.maximum_depth < 2 || options.maximum_depth > 16 ||
        options.threads < 1 || options.threads > 256 ||
        !(options.radius >= 1.0L) || !std::isfinite(options.radius)) {
        throw std::invalid_argument(
            "local-sld-response-diagonal requires states, a certificate, and valid max-depth/threads/radius");
    }
    return options;
}

void LocalSldResponseDiagonalCli::print_help(std::ostream& out) {
    out << "Local SLD cutoff-diagonal response options:\n"
        << "  --state PATH         add one cutoff state (repeat in cutoff order)\n"
        << "  --max-depth N        construct at most N response orders\n"
        << "  --radius R           weighted l1 radius r >= 1 (default 1.25)\n"
        << "  --threads N          direct bilinear-kernel workers\n"
        << "  --certificate PATH   write English JSON diagonal report\n";
}

int LocalSldResponseDiagonalCli::run(
    const LocalSldResponseDiagonalOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure("direct", options.threads);
    const SpectralDynamics dynamics(galerkin);
    std::vector<std::pair<std::string, SpectralState>> states;
    states.reserve(options.state_paths.size());
    for (const std::string& path : options.state_paths) {
        states.emplace_back(path, SpectralStateReader::read_tsv(path));
    }
    const LocalSldResponseDiagonalReport report =
        LocalSldResponseDiagonal::analyze(
            dynamics, states, options.maximum_depth, options.radius);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "response diagonal radius="
        << static_cast<double>(report.radius) << '\n'
        << "cutoff,safe_orders,weighted_l1,weighted_tail,last_safe,first_truncated,common_l2,quadratic_ratio\n";
    for (const LocalSldResponseDiagonalRow& row : report.rows) {
        out << row.cutoff << ',' << row.safe_order_count << ','
            << static_cast<double>(row.safe_weighted_l1) << ','
            << static_cast<double>(row.safe_weighted_tail_after_cubic)
            << ',' << static_cast<double>(row.last_safe_coefficient)
            << ',' << static_cast<double>(row.first_truncated_coefficient)
            << ',' << static_cast<double>(
                   row.common_coefficient_l2_difference)
            << ',' << static_cast<double>(row.quadratic_majorant_ratio)
            << '\n';
    }
    out << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
