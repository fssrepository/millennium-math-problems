#include "local_sld_response_family.hpp"

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
#include <utility>

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

void write_string_array(
    std::ostream& output,
    const std::vector<std::string>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        output << '"' << values[index] << '"';
        if (index + 1 != values.size()) {
            output << ", ";
        }
    }
    output << ']';
}

void write_certificate(
    const LocalSldResponseFamilyReport& report,
    const LocalSldResponseFamilyOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write response-family certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-response-family-v1\",\n"
        << "  \"depth\": " << report.depth << ",\n"
        << "  \"included_transverse_two_one_one\": "
        << (options.include_transverse_two_one_one ? "true" : "false")
        << ",\n"
        << "  \"included_three_one_zero_orbits\": "
        << (options.include_three_one_zero_orbits ? "true" : "false")
        << ",\n"
        << "  \"maximum_coefficient_l2_difference\": "
        << static_cast<double>(report.maximum_coefficient_l2_difference)
        << ",\n"
        << "  \"maximum_projection_residual\": "
        << static_cast<double>(report.maximum_projection_residual)
        << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const LocalSldResponseFamilyRow& row = report.rows[index];
        output << "    {\"state_path\": \"" << row.state_path
            << "\", \"cutoff\": " << row.cutoff
            << ", \"axis_response_energy\": "
            << static_cast<double>(row.axis_response_energy)
            << ", \"cubic_response_energy\": "
            << static_cast<double>(row.cubic_response_energy)
            << ", \"finite_response_energy\": "
            << static_cast<double>(row.finite_response_energy)
            << ", \"orbit_energy_increment\": "
            << static_cast<double>(row.orbit_energy_increment)
            << ", \"final_projection_energy\": "
            << static_cast<double>(row.final_projection_energy)
            << ", \"final_projection_residual\": "
            << static_cast<double>(row.final_projection_residual)
            << ", \"coefficient_l2_difference\": "
            << static_cast<double>(row.coefficient_l2_difference)
            << ", \"labels\": ";
        write_string_array(output, row.labels);
        output << ", \"coefficients\": ";
        write_real_array(output, row.coefficients);
        output << ", \"coefficient_energies\": ";
        write_real_array(output, row.coefficient_energies);
        output << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_cutoff_family_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldResponseFamilyReport LocalSldResponseFamily::analyze(
    const SpectralDynamics& dynamics,
    const std::vector<std::pair<std::string, SpectralState>>& states,
    int depth,
    bool include_transverse_two_one_one,
    bool include_three_one_zero_orbits) {
    if (states.empty()) {
        throw std::invalid_argument(
            "response family requires at least one state");
    }
    LocalSldResponseFamilyReport report;
    report.depth = depth;
    for (const auto& [path, state] : states) {
        const LocalSldResponseHierarchyReport hierarchy =
            LocalSldResponseHierarchy::analyze(
                dynamics, state, depth,
                include_transverse_two_one_one,
                include_three_one_zero_orbits);
        LocalSldResponseFamilyRow row;
        row.state_path = path;
        row.cutoff = hierarchy.cutoff;
        for (const LocalSldResponseHierarchyRow& entry : hierarchy.rows) {
            row.labels.push_back(entry.label);
            row.coefficients.push_back(entry.coefficient);
            row.coefficient_energies.push_back(entry.coefficient_energy);
            if (entry.order == 1) {
                row.axis_response_energy =
                    entry.cumulative_projection_energy;
            }
            if (entry.order == 2) {
                row.cubic_response_energy =
                    entry.cumulative_projection_energy;
            }
            if (entry.order == hierarchy.constructed_depth - 1) {
                row.finite_response_energy =
                    entry.cumulative_projection_energy;
            }
        }
        row.final_projection_energy =
            hierarchy.final_projection_energy;
        row.orbit_energy_increment = row.final_projection_energy -
            row.finite_response_energy;
        row.final_projection_residual =
            hierarchy.final_projection_residual;
        if (!report.rows.empty()) {
            const LocalSldResponseFamilyRow& previous =
                report.rows.back();
            const std::size_t common = std::min(
                previous.coefficients.size(), row.coefficients.size());
            SpectralReal difference2 = 0.0L;
            for (std::size_t index = 0; index < common; ++index) {
                const SpectralReal difference = row.coefficients[index] -
                    previous.coefficients[index];
                difference2 += difference * difference;
            }
            row.coefficient_l2_difference = std::sqrt(difference2);
        }
        report.maximum_coefficient_l2_difference = std::max(
            report.maximum_coefficient_l2_difference,
            row.coefficient_l2_difference);
        report.maximum_projection_residual = std::max(
            report.maximum_projection_residual,
            row.final_projection_residual);
        report.rows.push_back(std::move(row));
    }
    return report;
}

LocalSldResponseFamilyOptions LocalSldResponseFamilyCli::parse(
    int argc, char** argv, int first) {
    LocalSldResponseFamilyOptions options;
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
        } else if (name == "--depth") {
            options.depth = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--include-211-transverse") {
            options.include_transverse_two_one_one = true;
        } else if (name == "--include-310-orbits") {
            options.include_three_one_zero_orbits = true;
        } else {
            throw std::invalid_argument(
                "unknown local-sld-response-family option: " + name);
        }
    }
    if (options.state_paths.empty() || options.certificate_path.empty() ||
        options.depth < 2 || options.depth > 16 ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "local-sld-response-family requires repeated --state, --certificate, and valid depth/threads");
    }
    return options;
}

void LocalSldResponseFamilyCli::print_help(std::ostream& out) {
    out << "Local SLD response-family options:\n"
        << "  --state PATH         add one cutoff state (repeat in cutoff order)\n"
        << "  --depth N            response recursion depth (default 16)\n"
        << "  --threads N          direct bilinear-kernel workers\n"
        << "  --include-211-transverse add the second (2,1,1) polarization\n"
        << "  --include-310-orbits add both oriented (3,1,0) orbits\n"
        << "  --certificate PATH   write English JSON family report\n";
}

int LocalSldResponseFamilyCli::run(
    const LocalSldResponseFamilyOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure("direct", options.threads);
    const SpectralDynamics dynamics(galerkin);
    std::vector<std::pair<std::string, SpectralState>> states;
    states.reserve(options.state_paths.size());
    for (const std::string& path : options.state_paths) {
        states.emplace_back(path, SpectralStateReader::read_tsv(path));
    }
    const LocalSldResponseFamilyReport report =
        LocalSldResponseFamily::analyze(
            dynamics, states, options.depth,
            options.include_transverse_two_one_one,
            options.include_three_one_zero_orbits);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "response family depth=" << report.depth << '\n'
        << "cutoff,axis_response,cubic_response,finite_response,orbits,projection,residual,coefficient_difference\n";
    for (const LocalSldResponseFamilyRow& row : report.rows) {
        out << row.cutoff << ','
            << static_cast<double>(row.axis_response_energy) << ','
            << static_cast<double>(row.cubic_response_energy) << ','
            << static_cast<double>(row.finite_response_energy) << ','
            << static_cast<double>(row.orbit_energy_increment) << ','
            << static_cast<double>(row.final_projection_energy) << ','
            << static_cast<double>(row.final_projection_residual) << ','
            << static_cast<double>(row.coefficient_l2_difference) << '\n';
    }
    out << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
