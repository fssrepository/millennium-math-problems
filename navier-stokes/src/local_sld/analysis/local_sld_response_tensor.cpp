#include "local_sld_response_tensor.hpp"

#include "local_sld_cyclic_basis.hpp"
#include "spectral_galerkin.hpp"

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

int highest_active_shell(const SpectralState& state) {
    SpectralReal maximum_mode_energy = 0.0L;
    for (const ComplexVector& value : state.velocity) {
        maximum_mode_energy = std::max(
            maximum_mode_energy,
            std::real(dot_hermitian(value, value)));
    }
    const SpectralReal threshold = 1e-24L * maximum_mode_energy;
    int shell = 0;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        if (std::real(dot_hermitian(
                state.velocity[mode], state.velocity[mode])) <= threshold) {
            continue;
        }
        const WaveVector wave = state.waves[mode];
        shell = std::max(shell, std::max(
            {std::abs(wave.x), std::abs(wave.y), std::abs(wave.z)}));
    }
    return shell;
}

SpectralReal maximum_gram_error(
    const std::vector<SpectralState>& basis) {
    SpectralReal error = 0.0L;
    for (std::size_t left = 0; left < basis.size(); ++left) {
        for (std::size_t right = 0; right < basis.size(); ++right) {
            const SpectralReal expected = left == right ? 1.0L : 0.0L;
            error = std::max(error, std::abs(
                LocalSldCyclicBasis::pairing(
                    basis[left].velocity, basis[right].velocity) -
                expected));
        }
    }
    return error;
}

void write_certificate(
    const LocalSldResponseTensorReport& report,
    const LocalSldResponseTensorOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write response-tensor certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-response-tensor-v2\",\n"
        << "  \"definition\": \"exact direct-triad coefficients <b_m,B(b_i,b_j)> on the cutoff-diagonal response basis\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"depth\": " << report.depth << ",\n"
        << "  \"input_radius\": "
        << static_cast<double>(report.input_radius) << ",\n"
        << "  \"output_radius\": "
        << static_cast<double>(report.output_radius) << ",\n"
        << "  \"entry_tolerance\": "
        << static_cast<double>(report.entry_tolerance) << ",\n"
        << "  \"boundary_free_depth\": "
        << (report.boundary_free_depth ? "true" : "false") << ",\n"
        << "  \"maximum_gram_error\": "
        << static_cast<double>(report.maximum_gram_error) << ",\n"
        << "  \"maximum_projected_bilinear_constant\": "
        << static_cast<double>(
               report.maximum_projected_bilinear_constant)
        << ",\n"
        << "  \"maximum_complement_norm\": "
        << static_cast<double>(report.maximum_complement_norm) << ",\n"
        << "  \"maximum_complement_fraction\": "
        << static_cast<double>(report.maximum_complement_fraction) << ",\n"
        << "  \"maximum_input_weighted_complement_norm\": "
        << static_cast<double>(
               report.maximum_input_weighted_complement_norm)
        << ",\n"
        << "  \"maximum_norm_reconstruction_error\": "
        << static_cast<double>(report.maximum_norm_reconstruction_error)
        << ",\n"
        << "  \"retained_tensor_entries\": "
        << report.retained_tensor_entries << ",\n"
        << "  \"pairs\": [\n";
    for (std::size_t pair_index = 0;
         pair_index < report.pairs.size(); ++pair_index) {
        const LocalSldResponseTensorPair& pair =
            report.pairs[pair_index];
        output << "    {\"left_order\": " << pair.left_order
            << ", \"right_order\": " << pair.right_order
            << ", \"left_highest_shell\": "
            << pair.left_highest_shell
            << ", \"right_highest_shell\": "
            << pair.right_highest_shell
            << ", \"full_norm\": "
            << static_cast<double>(pair.full_norm)
            << ", \"projected_norm\": "
            << static_cast<double>(pair.projected_norm)
            << ", \"complement_norm\": "
            << static_cast<double>(pair.complement_norm)
            << ", \"complement_fraction\": "
            << static_cast<double>(pair.complement_fraction)
            << ", \"weighted_projected_ratio\": "
            << static_cast<double>(pair.weighted_projected_ratio)
            << ", \"input_weighted_complement_norm\": "
            << static_cast<double>(
                   pair.input_weighted_complement_norm)
            << ", \"entries\": [";
        for (std::size_t entry_index = 0;
             entry_index < pair.entries.size(); ++entry_index) {
            const LocalSldResponseTensorEntry& entry =
                pair.entries[entry_index];
            output << "{\"output_order\": " << entry.output_order
                << ", \"coefficient\": "
                << static_cast<double>(entry.coefficient) << '}';
            if (entry_index + 1 != pair.entries.size()) {
                output << ", ";
            }
        }
        output << "]}"
            << (pair_index + 1 == report.pairs.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_tensor_is_not_a_proof\": true,\n"
        << "  \"remaining_requirement\": \"prove cutoff-uniform projected and complement bilinear bounds in the weighted response norm\"\n"
        << "}\n";
}

}  // namespace

LocalSldResponseTensorReport LocalSldResponseTensor::analyze(
    const SpectralDynamics& dynamics,
    int cutoff,
    int depth,
    SpectralReal input_radius,
    SpectralReal output_radius,
    SpectralReal entry_tolerance) {
    if (cutoff < 2 || cutoff > 12 || depth < 2 || depth > 16 ||
        depth > cutoff + 1 || !(input_radius >= 1.0L) ||
        !(output_radius >= 1.0L) ||
        !(input_radius >= output_radius) ||
        !std::isfinite(input_radius) ||
        !std::isfinite(output_radius) ||
        !(entry_tolerance >= 0.0L) ||
        !std::isfinite(entry_tolerance)) {
        throw std::invalid_argument(
            "response tensor requires cutoff 2..12, depth 2..min(16,K+1), radius >= 1, and finite tolerance");
    }
    LocalSldResponseTensorReport report;
    report.cutoff = cutoff;
    report.depth = depth;
    report.input_radius = input_radius;
    report.output_radius = output_radius;
    report.entry_tolerance = entry_tolerance;
    report.boundary_free_depth = depth <= cutoff + 1;
    const std::vector<SpectralState> basis =
        LocalSldResponseHierarchy::build(dynamics, cutoff, depth);
    if (static_cast<int>(basis.size()) != depth) {
        throw std::runtime_error(
            "response recursion terminated before requested tensor depth");
    }
    report.maximum_gram_error = maximum_gram_error(basis);
    std::vector<int> shells;
    shells.reserve(basis.size());
    for (const SpectralState& state : basis) {
        shells.push_back(highest_active_shell(state));
    }
    report.pairs.reserve(basis.size() * basis.size());
    for (int left = 0; left < depth; ++left) {
        for (int right = 0; right < depth; ++right) {
            const SpectralIncrement product =
                dynamics.advection_bilinear_direct_partition(
                    basis.front(),
                    basis[static_cast<std::size_t>(left)].velocity,
                    basis[static_cast<std::size_t>(right)].velocity,
                    TriadPartition::all);
            LocalSldResponseTensorPair pair;
            pair.left_order = left;
            pair.right_order = right;
            pair.left_highest_shell = shells[
                static_cast<std::size_t>(left)];
            pair.right_highest_shell = shells[
                static_cast<std::size_t>(right)];
            const SpectralReal full_norm2 = std::max(
                0.0L,
                LocalSldCyclicBasis::pairing(product, product));
            SpectralReal projected_norm2 = 0.0L;
            SpectralReal weighted_sum = 0.0L;
            for (int output_order = 0;
                 output_order < depth; ++output_order) {
                const SpectralReal coefficient =
                    LocalSldCyclicBasis::pairing(
                        basis[static_cast<std::size_t>(output_order)]
                            .velocity,
                        product);
                projected_norm2 += coefficient * coefficient;
                weighted_sum += std::pow(
                    output_radius, output_order) *
                    std::abs(coefficient);
                if (std::abs(coefficient) > entry_tolerance) {
                    pair.entries.push_back(
                        LocalSldResponseTensorEntry{
                            output_order, coefficient});
                    ++report.retained_tensor_entries;
                }
            }
            const SpectralReal complement_norm2 = std::max(
                0.0L, full_norm2 - projected_norm2);
            pair.full_norm = std::sqrt(full_norm2);
            pair.projected_norm = std::sqrt(
                std::max(0.0L, projected_norm2));
            pair.complement_norm = std::sqrt(complement_norm2);
            pair.complement_fraction = pair.full_norm > 0.0L
                ? pair.complement_norm / pair.full_norm
                : 0.0L;
            pair.weighted_projected_ratio = weighted_sum /
                std::pow(input_radius, left + right);
            pair.input_weighted_complement_norm =
                pair.complement_norm /
                std::pow(input_radius, left + right);
            const SpectralReal reconstruction_error = std::abs(
                full_norm2 - projected_norm2 - complement_norm2) /
                std::max(full_norm2, 1e-30L);
            report.maximum_norm_reconstruction_error = std::max(
                report.maximum_norm_reconstruction_error,
                reconstruction_error);
            report.maximum_projected_bilinear_constant = std::max(
                report.maximum_projected_bilinear_constant,
                pair.weighted_projected_ratio);
            report.maximum_complement_norm = std::max(
                report.maximum_complement_norm,
                pair.complement_norm);
            report.maximum_complement_fraction = std::max(
                report.maximum_complement_fraction,
                pair.complement_fraction);
            report.maximum_input_weighted_complement_norm = std::max(
                report.maximum_input_weighted_complement_norm,
                pair.input_weighted_complement_norm);
            report.pairs.push_back(std::move(pair));
        }
    }
    return report;
}

LocalSldResponseTensorOptions LocalSldResponseTensorCli::parse(
    int argc, char** argv, int first) {
    LocalSldResponseTensorOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--cutoff") {
            options.cutoff = std::stoi(next(index, name));
        } else if (name == "--depth") {
            options.depth = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--radius") {
            const SpectralReal radius = std::stold(next(index, name));
            options.input_radius = radius;
            options.output_radius = radius;
        } else if (name == "--input-radius") {
            options.input_radius = std::stold(next(index, name));
        } else if (name == "--output-radius") {
            options.output_radius = std::stold(next(index, name));
        } else if (name == "--tolerance") {
            options.entry_tolerance = std::stold(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-sld-response-tensor option: " + name);
        }
    }
    if (options.certificate_path.empty() || options.cutoff < 2 ||
        options.cutoff > 12 || options.depth < 2 ||
        options.depth > 16 || options.depth > options.cutoff + 1 ||
        options.threads < 1 || options.threads > 256 ||
        !(options.input_radius >= options.output_radius) ||
        !(options.output_radius >= 1.0L) ||
        !std::isfinite(options.input_radius) ||
        !std::isfinite(options.output_radius) ||
        !(options.entry_tolerance >= 0.0L) ||
        !std::isfinite(options.entry_tolerance)) {
        throw std::invalid_argument(
            "local-sld-response-tensor requires a certificate and valid cutoff/depth/threads/radius/tolerance");
    }
    return options;
}

void LocalSldResponseTensorCli::print_help(std::ostream& out) {
    out << "Local SLD response interaction-tensor options:\n"
        << "  --cutoff K           Galerkin cutoff in 2..12\n"
        << "  --depth N            boundary-free basis depth <= K+1\n"
        << "  --input-radius R     stronger input radius R >= r\n"
        << "  --output-radius r    weaker output radius r >= 1\n"
        << "  --radius r           shorthand setting both radii equally\n"
        << "  --tolerance X        sparse-entry output threshold\n"
        << "  --threads N          direct bilinear-kernel workers\n"
        << "  --certificate PATH   write English JSON tensor ledger\n";
}

int LocalSldResponseTensorCli::run(
    const LocalSldResponseTensorOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure("direct", options.threads);
    const SpectralDynamics dynamics(galerkin);
    const LocalSldResponseTensorReport report =
        LocalSldResponseTensor::analyze(
            dynamics, options.cutoff, options.depth,
            options.input_radius, options.output_radius,
            options.entry_tolerance);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "response tensor cutoff=" << report.cutoff
        << " depth=" << report.depth
        << " radii=" << static_cast<double>(report.input_radius)
        << "->" << static_cast<double>(report.output_radius)
        << " entries=" << report.retained_tensor_entries
        << " projected_constant="
        << static_cast<double>(
               report.maximum_projected_bilinear_constant)
        << " complement_norm="
        << static_cast<double>(report.maximum_complement_norm)
        << " complement_fraction="
        << static_cast<double>(report.maximum_complement_fraction)
        << " weighted_complement="
        << static_cast<double>(
               report.maximum_input_weighted_complement_norm)
        << " Gram_error="
        << static_cast<double>(report.maximum_gram_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
