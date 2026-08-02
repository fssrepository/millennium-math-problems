#include "local_sld_response_hierarchy.hpp"

#include "local_sld_cyclic_basis.hpp"
#include "local_sld_cyclic_orbit_basis.hpp"
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

void add_increment(
    SpectralIncrement& target,
    const SpectralIncrement& source) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "response hierarchy increment layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += source[mode][component];
        }
    }
}

int highest_active_shell(const SpectralState& state) {
    SpectralReal maximum_mode_energy = 0.0L;
    for (const ComplexVector& value : state.velocity) {
        maximum_mode_energy = std::max(
            maximum_mode_energy,
            std::real(dot_hermitian(value, value)));
    }
    const SpectralReal threshold = 1e-24L * maximum_mode_energy;
    int result = 0;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        if (std::real(dot_hermitian(
                state.velocity[mode], state.velocity[mode])) <= threshold) {
            continue;
        }
        const WaveVector wave = state.waves[mode];
        result = std::max(
            result,
            std::max({std::abs(wave.x), std::abs(wave.y),
                      std::abs(wave.z)}));
    }
    return result;
}

SpectralReal maximum_gram_error(
    const std::vector<SpectralState>& basis) {
    SpectralReal result = 0.0L;
    for (std::size_t left = 0; left < basis.size(); ++left) {
        for (std::size_t right = 0; right < basis.size(); ++right) {
            const SpectralReal expected = left == right ? 1.0L : 0.0L;
            result = std::max(result, std::abs(
                LocalSldCyclicBasis::pairing(
                    basis[left].velocity, basis[right].velocity) -
                expected));
        }
    }
    return result;
}

void append_orthogonalized(
    const SpectralDynamics& dynamics,
    std::vector<SpectralState>& basis,
    SpectralState candidate) {
    for (const SpectralState& previous : basis) {
        const SpectralReal weight = LocalSldCyclicBasis::pairing(
            previous.velocity, candidate.velocity);
        for (std::size_t mode = 0;
             mode < candidate.velocity.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                candidate.velocity[mode][component] -=
                    weight * previous.velocity[mode][component];
            }
        }
    }
    dynamics.enforce_constraints(candidate);
    SpectralStateOps::normalize_energy(candidate);
    basis.push_back(std::move(candidate));
}

void write_certificate(
    const LocalSldResponseHierarchyReport& report,
    const LocalSldResponseHierarchyOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write response hierarchy certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-response-hierarchy-v1\",\n"
        << "  \"construction\": \"Gram-Schmidt orthogonalized coefficients of the quadratic response recursion sum_{i+j=n-1} B(b_i,b_j)\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"residual_state_path\": \""
        << options.residual_state_path << "\",\n"
        << "  \"projected_state_path\": \""
        << options.projected_state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"requested_depth\": "
        << report.requested_depth << ",\n"
        << "  \"constructed_depth\": "
        << report.constructed_depth << ",\n"
        << "  \"included_transverse_two_one_one\": "
        << (report.included_transverse_two_one_one ? "true" : "false")
        << ",\n"
        << "  \"included_three_one_zero_orbits\": "
        << (report.included_three_one_zero_orbits ? "true" : "false")
        << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"reference_energy\": "
        << static_cast<double>(report.reference_energy) << ",\n"
        << "  \"maximum_gram_error\": "
        << static_cast<double>(report.maximum_gram_error) << ",\n"
        << "  \"final_projection_energy\": "
        << static_cast<double>(report.final_projection_energy) << ",\n"
        << "  \"final_projection_residual\": "
        << static_cast<double>(report.final_projection_residual) << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const LocalSldResponseHierarchyRow& row = report.rows[index];
        output << "    {\"order\": " << row.order
            << ", \"label\": \"" << row.label << "\""
            << ", \"coefficient\": "
            << static_cast<double>(row.coefficient)
            << ", \"coefficient_energy\": "
            << static_cast<double>(row.coefficient_energy)
            << ", \"cumulative_projection_energy\": "
            << static_cast<double>(row.cumulative_projection_energy)
            << ", \"projection_residual\": "
            << static_cast<double>(row.projection_residual)
            << ", \"highest_active_shell\": "
            << row.highest_active_shell << "}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_projection_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

std::vector<SpectralState> LocalSldResponseHierarchy::build(
    const SpectralDynamics& dynamics,
    int cutoff,
    int depth) {
    if (cutoff < 2 || cutoff > 12 || depth < 2 || depth > 16) {
        throw std::invalid_argument(
            "response hierarchy requires cutoff 2..12 and depth 2..16");
    }
    std::vector<SpectralState> basis;
    basis.reserve(static_cast<std::size_t>(depth));
    basis.push_back(LocalSldCyclicBasis::axis_state(cutoff));
    basis.push_back(LocalSldCyclicBasis::response_state(
        dynamics, basis.front()));
    for (int order = 2; order < depth; ++order) {
        SpectralState candidate = basis.front();
        for (ComplexVector& value : candidate.velocity) {
            value = {};
        }
        for (int left = 0; left < order; ++left) {
            const int right = order - 1 - left;
            add_increment(
                candidate.velocity,
                dynamics.advection_bilinear_direct_partition(
                    candidate,
                    basis[static_cast<std::size_t>(left)].velocity,
                    basis[static_cast<std::size_t>(right)].velocity,
                    TriadPartition::all));
        }
        dynamics.enforce_constraints(candidate);
        for (const SpectralState& previous : basis) {
            const SpectralReal weight = LocalSldCyclicBasis::pairing(
                previous.velocity, candidate.velocity);
            for (std::size_t mode = 0;
                 mode < candidate.velocity.size(); ++mode) {
                for (std::size_t component = 0; component < 3;
                     ++component) {
                    candidate.velocity[mode][component] -=
                        weight * previous.velocity[mode][component];
                }
            }
        }
        dynamics.enforce_constraints(candidate);
        if (SpectralStateOps::energy(candidate) < 1e-28L) {
            break;
        }
        SpectralStateOps::normalize_energy(candidate);
        basis.push_back(std::move(candidate));
    }
    return basis;
}

LocalSldResponseHierarchyReport LocalSldResponseHierarchy::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& reference,
    int depth,
    bool include_transverse_two_one_one,
    bool include_three_one_zero_orbits) {
    LocalSldResponseHierarchyReport report;
    report.cutoff = SpectralStateOps::cutoff(reference);
    report.requested_depth = depth;
    report.reference_energy = SpectralStateOps::energy(reference);
    std::vector<SpectralState> basis = build(
        dynamics, report.cutoff, depth);
    report.constructed_depth = static_cast<int>(basis.size());
    report.included_transverse_two_one_one =
        include_transverse_two_one_one;
    report.included_three_one_zero_orbits =
        include_three_one_zero_orbits;
    if (include_transverse_two_one_one) {
        append_orthogonalized(
            dynamics, basis,
            LocalSldCyclicOrbitBasis::transverse_two_one_one(
                report.cutoff));
    }
    if (include_three_one_zero_orbits) {
        append_orthogonalized(
            dynamics, basis,
            LocalSldCyclicOrbitBasis::forward_three_one_zero(
                report.cutoff));
        append_orthogonalized(
            dynamics, basis,
            LocalSldCyclicOrbitBasis::backward_three_one_zero(
                report.cutoff));
    }
    report.maximum_gram_error = maximum_gram_error(basis);
    report.residual_state = reference;
    SpectralReal cumulative = 0.0L;
    for (std::size_t order = 0; order < basis.size(); ++order) {
        LocalSldResponseHierarchyRow row;
        row.order = static_cast<int>(order);
        if (order < static_cast<std::size_t>(report.constructed_depth)) {
            row.label = "response-order-" + std::to_string(order);
        } else {
            std::size_t extra = order - static_cast<std::size_t>(
                report.constructed_depth);
            if (include_transverse_two_one_one) {
                if (extra == 0) {
                    row.label = "transverse-(2,1,1)-orbit";
                } else {
                    --extra;
                }
            }
            if (row.label.empty()) {
                row.label = extra == 0
                    ? "forward-(3,1,0)-orbit"
                    : "backward-(3,1,0)-orbit";
            }
        }
        row.coefficient = LocalSldCyclicBasis::pairing(
            reference.velocity, basis[order].velocity);
        row.coefficient_energy =
            row.coefficient * row.coefficient / report.reference_energy;
        cumulative += row.coefficient_energy;
        row.cumulative_projection_energy = cumulative;
        row.projection_residual = std::sqrt(std::max(
            0.0L, 1.0L - cumulative));
        row.highest_active_shell = highest_active_shell(basis[order]);
        report.rows.push_back(row);
        for (std::size_t mode = 0;
             mode < report.residual_state.velocity.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                report.residual_state.velocity[mode][component] -=
                    row.coefficient *
                    basis[order].velocity[mode][component];
            }
        }
    }
    report.final_projection_energy = cumulative;
    report.final_projection_residual = std::sqrt(std::max(
        0.0L, 1.0L - cumulative));
    dynamics.enforce_constraints(report.residual_state);
    report.projected_state = reference;
    for (std::size_t mode = 0;
         mode < report.projected_state.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            report.projected_state.velocity[mode][component] -=
                report.residual_state.velocity[mode][component];
        }
    }
    dynamics.enforce_constraints(report.projected_state);
    return report;
}

LocalSldResponseHierarchyOptions LocalSldResponseHierarchyCli::parse(
    int argc, char** argv, int first) {
    LocalSldResponseHierarchyOptions options;
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
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--residual-state") {
            options.residual_state_path = next(index, name);
        } else if (name == "--projected-state") {
            options.projected_state_path = next(index, name);
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
                "unknown local-sld-response-hierarchy option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.depth < 2 || options.depth > 16 ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "local-sld-response-hierarchy requires --state, --certificate, and valid depth/threads");
    }
    return options;
}

void LocalSldResponseHierarchyCli::print_help(std::ostream& out) {
    out << "Local SLD quadratic response hierarchy options:\n"
        << "  --state PATH         reference Fourier TSV\n"
        << "  --depth N            response orders to construct (default 6)\n"
        << "  --threads N          direct bilinear-kernel workers\n"
        << "  --include-211-transverse add the second cyclic polarization orbit\n"
        << "  --include-310-orbits add both oriented cyclic (3,1,0) orbits\n"
        << "  --certificate PATH   write English JSON projection report\n";
    out << "  --residual-state PATH write the unnormalized projection residual\n";
    out << "  --projected-state PATH write the unnormalized projected state\n";
}

int LocalSldResponseHierarchyCli::run(
    const LocalSldResponseHierarchyOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure("direct", options.threads);
    const SpectralDynamics dynamics(galerkin);
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    const LocalSldResponseHierarchyReport report =
        LocalSldResponseHierarchy::analyze(
            dynamics, state, options.depth,
            options.include_transverse_two_one_one,
            options.include_three_one_zero_orbits);
    write_certificate(report, options);
    if (!options.residual_state_path.empty()) {
        SpectralStateWriter::write_tsv(
            options.residual_state_path, report.residual_state,
            "quadratic response-hierarchy projection residual; candidate_lemma_proved=false");
    }
    if (!options.projected_state_path.empty()) {
        SpectralStateWriter::write_tsv(
            options.projected_state_path, report.projected_state,
            "quadratic response hierarchy plus cyclic orbit projection; candidate_lemma_proved=false");
    }
    out << std::setprecision(12)
        << "response hierarchy cutoff=" << report.cutoff
        << " depth=" << report.constructed_depth
        << " Gram_error="
        << static_cast<double>(report.maximum_gram_error)
        << " projection_energy="
        << static_cast<double>(report.final_projection_energy)
        << " residual="
        << static_cast<double>(report.final_projection_residual) << '\n'
        << "order,label,coefficient,coefficient_energy,cumulative,residual,shell\n";
    for (const LocalSldResponseHierarchyRow& row : report.rows) {
        out << row.order << ',' << row.label << ','
            << static_cast<double>(row.coefficient) << ','
            << static_cast<double>(row.coefficient_energy) << ','
            << static_cast<double>(row.cumulative_projection_energy) << ','
            << static_cast<double>(row.projection_residual) << ','
            << row.highest_active_shell << '\n';
    }
    out << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
