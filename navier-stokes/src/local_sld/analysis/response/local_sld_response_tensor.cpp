#include "local_sld_response_tensor.hpp"

#include "local_sld_cyclic_basis.hpp"
#include "local_sld_response_basis.hpp"
#include "spectral_galerkin.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace lemma {
namespace {

SpectralState quadratic_product_state(
    const SpectralDynamics& dynamics,
    const std::vector<LocalSldResponseBasisElement>& basis,
    std::size_t left,
    std::size_t right) {
    SpectralState candidate = basis.front().state;
    candidate.velocity = dynamics.advection_bilinear_direct_partition(
        candidate,
        basis[left].state.velocity,
        basis[right].state.velocity,
        TriadPartition::all);
    dynamics.enforce_constraints(candidate);
    return candidate;
}

SpectralState orthogonal_complement_candidate(
    const SpectralDynamics& dynamics,
    const std::vector<LocalSldResponseBasisElement>& basis,
    SpectralState candidate) {
    for (const LocalSldResponseBasisElement& output : basis) {
        const SpectralReal coefficient = LocalSldCyclicBasis::pairing(
            output.state.velocity, candidate.velocity);
        for (std::size_t mode = 0;
             mode < candidate.velocity.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                candidate.velocity[mode][component] -=
                    coefficient * output.state.velocity[mode][component];
            }
        }
    }
    dynamics.enforce_constraints(candidate);
    return candidate;
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
        << "  \"schema\": \"navier-stokes-local-sld-response-tensor-v8\",\n"
        << "  \"definition\": \"exact direct-triad coefficients <b_m,B(b_i,b_j)> on the cutoff-diagonal response basis\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"depth\": " << report.depth << ",\n"
        << "  \"input_radius\": "
        << static_cast<double>(report.input_radius) << ",\n"
        << "  \"output_radius\": "
        << static_cast<double>(report.output_radius) << ",\n"
        << "  \"entry_tolerance\": "
        << static_cast<double>(report.entry_tolerance) << ",\n"
        << "  \"included_transverse_two_one_one\": "
        << (report.included_transverse_two_one_one ? "true" : "false")
        << ",\n"
        << "  \"included_three_one_zero_orbits\": "
        << (report.included_three_one_zero_orbits ? "true" : "false")
        << ",\n"
        << "  \"basis_size\": " << report.basis_size << ",\n"
        << "  \"closure_extensions_requested\": "
        << report.closure_extensions_requested << ",\n"
        << "  \"closure_extensions_constructed\": "
        << report.closure_extensions_constructed << ",\n"
        << "  \"boundary_free_depth\": "
        << (report.boundary_free_depth ? "true" : "false") << ",\n"
        << "  \"maximum_gram_error\": "
        << static_cast<double>(report.maximum_gram_error) << ",\n"
        << "  \"maximum_projected_bilinear_constant\": "
        << static_cast<double>(
               report.maximum_projected_bilinear_constant)
        << ",\n"
        << "  \"maximum_projected_degree_block_constant\": "
        << static_cast<double>(
               report.maximum_projected_degree_block_constant)
        << ",\n"
        << "  \"axis_pair_candidate_constant\": "
        << static_cast<double>(report.axis_pair_candidate_constant)
        << ",\n"
        << "  \"projected_constant_over_axis_candidate\": "
        << static_cast<double>(
               report.projected_constant_over_axis_candidate)
        << ",\n"
        << "  \"axis_pair_candidate_survives\": "
        << (report.axis_pair_candidate_survives ? "true" : "false")
        << ",\n"
        << "  \"maximum_complement_norm\": "
        << static_cast<double>(report.maximum_complement_norm) << ",\n"
        << "  \"maximum_complement_fraction\": "
        << static_cast<double>(report.maximum_complement_fraction) << ",\n"
        << "  \"maximum_input_weighted_complement_norm\": "
        << static_cast<double>(
               report.maximum_input_weighted_complement_norm)
        << ",\n"
        << "  \"maximum_shell_weighted_complement_ratio\": "
        << static_cast<double>(
               report.maximum_shell_weighted_complement_ratio)
        << ",\n"
        << "  \"maximum_input_weighted_h1_complement_norm\": "
        << static_cast<double>(
               report.maximum_input_weighted_h1_complement_norm)
        << ",\n"
        << "  \"maximum_input_weighted_h2_complement_norm\": "
        << static_cast<double>(
               report.maximum_input_weighted_h2_complement_norm)
        << ",\n"
        << "  \"finite_combined_weighted_bound\": "
        << static_cast<double>(report.finite_combined_weighted_bound)
        << ",\n"
        << "  \"maximum_norm_reconstruction_error\": "
        << static_cast<double>(report.maximum_norm_reconstruction_error)
        << ",\n"
        << "  \"maximum_shell_norm_reconstruction_error\": "
        << static_cast<double>(
               report.maximum_shell_norm_reconstruction_error)
        << ",\n"
        << "  \"maximum_output_degree_excess\": "
        << report.maximum_output_degree_excess << ",\n"
        << "  \"maximum_degree_excess_coefficient\": "
        << static_cast<double>(
               report.maximum_degree_excess_coefficient)
        << ",\n"
        << "  \"graded_support_closed\": "
        << (report.graded_support_closed ? "true" : "false")
        << ",\n"
        << "  \"retained_tensor_entries\": "
        << report.retained_tensor_entries << ",\n"
        << "  \"basis\": [\n";
    for (std::size_t index = 0; index < report.basis.size(); ++index) {
        const LocalSldResponseTensorBasisRow& row = report.basis[index];
        output << "    {\"basis_index\": " << row.basis_index
            << ", \"response_order\": " << row.response_order
            << ", \"highest_active_shell\": "
            << row.highest_active_shell
            << ", \"analytic_degree\": " << row.analytic_degree
            << ", \"label\": \"" << row.label
            << "\", \"scalar_response\": "
            << (row.scalar_response ? "true" : "false") << '}'
            << (index + 1 == report.basis.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"degree_offset_envelope\": [\n";
    for (std::size_t index = 0;
         index < report.degree_offset_envelope.size(); ++index) {
        const LocalSldResponseTensorDegreeOffsetRow& row =
            report.degree_offset_envelope[index];
        output << "    {\"degree_offset\": " << row.degree_offset
            << ", \"retained_entries\": " << row.retained_entries
            << ", \"maximum_absolute_coefficient\": "
            << static_cast<double>(row.maximum_absolute_coefficient)
            << ", \"maximum_weighted_entry\": "
            << static_cast<double>(row.maximum_weighted_entry)
            << ", \"maximum_pair_block_contribution\": "
            << static_cast<double>(
                   row.maximum_pair_block_contribution) << '}'
            << (index + 1 == report.degree_offset_envelope.size()
                ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"pairs\": [\n";
    for (std::size_t pair_index = 0;
         pair_index < report.pairs.size(); ++pair_index) {
        const LocalSldResponseTensorPair& pair =
            report.pairs[pair_index];
        output << "    {\"left_basis_index\": "
            << pair.left_basis_index
            << ", \"right_basis_index\": "
            << pair.right_basis_index
            << ", \"left_order\": " << pair.left_order
            << ", \"right_order\": " << pair.right_order
            << ", \"left_analytic_degree\": "
            << pair.left_analytic_degree
            << ", \"right_analytic_degree\": "
            << pair.right_analytic_degree
            << ", \"left_label\": \"" << pair.left_label
            << "\", \"right_label\": \"" << pair.right_label << '"'
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
            << ", \"weighted_projected_degree_block_ratio\": "
            << static_cast<double>(
                   pair.weighted_projected_degree_block_ratio)
            << ", \"input_weighted_complement_norm\": "
            << static_cast<double>(
                   pair.input_weighted_complement_norm)
            << ", \"shell_weighted_complement_ratio\": "
            << static_cast<double>(
                   pair.shell_weighted_complement_ratio)
            << ", \"h1_complement_norm\": "
            << static_cast<double>(pair.h1_complement_norm)
            << ", \"h2_complement_norm\": "
            << static_cast<double>(pair.h2_complement_norm)
            << ", \"input_weighted_h1_complement_norm\": "
            << static_cast<double>(
                   pair.input_weighted_h1_complement_norm)
            << ", \"input_weighted_h2_complement_norm\": "
            << static_cast<double>(
                   pair.input_weighted_h2_complement_norm)
            << ", \"shell_norm_reconstruction_error\": "
            << static_cast<double>(
                   pair.shell_norm_reconstruction_error)
            << ", \"maximum_output_degree_excess\": "
            << pair.maximum_output_degree_excess
            << ", \"maximum_degree_excess_coefficient\": "
            << static_cast<double>(
                   pair.maximum_degree_excess_coefficient)
            << ", \"entries\": [";
        for (std::size_t entry_index = 0;
             entry_index < pair.entries.size(); ++entry_index) {
            const LocalSldResponseTensorEntry& entry =
                pair.entries[entry_index];
            output << "{\"output_basis_index\": "
                << entry.output_basis_index
                << ", \"output_order\": " << entry.output_order
                << ", \"output_analytic_degree\": "
                << entry.output_analytic_degree
                << ", \"output_label\": \"" << entry.output_label
                << '"'
                << ", \"coefficient\": "
                << static_cast<double>(entry.coefficient) << '}';
            if (entry_index + 1 != pair.entries.size()) {
                output << ", ";
            }
        }
        output << "], \"complement_shells\": [";
        for (std::size_t shell_index = 0;
             shell_index < pair.complement_shells.size(); ++shell_index) {
            const LocalSldResponseTensorShell& shell =
                pair.complement_shells[shell_index];
            output << "{\"shell\": " << shell.shell
                << ", \"norm\": " << static_cast<double>(shell.norm)
                << ", \"complement_fraction\": "
                << static_cast<double>(shell.complement_fraction)
                << ", \"output_weighted_norm\": "
                << static_cast<double>(shell.output_weighted_norm)
                << '}';
            if (shell_index + 1 != pair.complement_shells.size()) {
                output << ", ";
            }
        }
        output << "]}"
            << (pair_index + 1 == report.pairs.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_tensor_is_not_a_proof\": true,\n"
        << "  \"remaining_requirement\": \"prove cutoff-uniform projected and shell-weighted transverse complement bounds compatible with the SLD derivative norm\"\n"
        << "}\n";
}

}  // namespace

LocalSldResponseTensorReport LocalSldResponseTensor::analyze(
    const SpectralDynamics& dynamics,
    int cutoff,
    int depth,
    SpectralReal input_radius,
    SpectralReal output_radius,
    SpectralReal entry_tolerance,
    bool include_transverse_two_one_one,
    bool include_three_one_zero_orbits,
    int closure_extensions) {
    if (cutoff < 2 || cutoff > 12 || depth < 2 || depth > 16 ||
        depth > cutoff + 1 || !(input_radius >= 1.0L) ||
        !(output_radius >= 1.0L) ||
        !(input_radius >= output_radius) ||
        !std::isfinite(input_radius) ||
        !std::isfinite(output_radius) ||
        !(entry_tolerance >= 0.0L) ||
        !std::isfinite(entry_tolerance) || closure_extensions < 0 ||
        closure_extensions > 16) {
        throw std::invalid_argument(
            "response tensor requires cutoff 2..12, depth 2..min(16,K+1), radius >= 1, and finite tolerance");
    }
    LocalSldResponseTensorReport report;
    report.cutoff = cutoff;
    report.depth = depth;
    report.input_radius = input_radius;
    report.output_radius = output_radius;
    report.entry_tolerance = entry_tolerance;
    report.included_transverse_two_one_one =
        include_transverse_two_one_one;
    report.included_three_one_zero_orbits =
        include_three_one_zero_orbits;
    report.closure_extensions_requested = closure_extensions;
    report.boundary_free_depth = depth <= cutoff + 1;
    std::vector<LocalSldResponseBasisElement> basis =
        LocalSldResponseBasis::build(
            dynamics, cutoff, depth,
            include_transverse_two_one_one,
            include_three_one_zero_orbits);
    const std::size_t expected_size = static_cast<std::size_t>(
        depth + (include_transverse_two_one_one ? 1 : 0) +
        (include_three_one_zero_orbits ? 2 : 0));
    if (basis.size() != expected_size) {
        throw std::runtime_error(
            "response recursion terminated before requested tensor depth");
    }
    for (int extension = 0;
         extension < closure_extensions; ++extension) {
        SpectralReal best_score = 0.0L;
        int best_degree = 0;
        std::size_t best_left = 0;
        std::size_t best_right = 0;
        SpectralState best_raw_candidate;
        for (std::size_t left = 0; left < basis.size(); ++left) {
            for (std::size_t right = 0; right < basis.size(); ++right) {
                SpectralState raw_candidate = quadratic_product_state(
                    dynamics, basis, left, right);
                SpectralState candidate = orthogonal_complement_candidate(
                    dynamics, basis, raw_candidate);
                const SpectralReal norm = std::sqrt(std::max(
                    0.0L, SpectralStateOps::energy(candidate)));
                const int degree = basis[left].analytic_degree +
                    basis[right].analytic_degree + 1;
                const SpectralReal score =
                    std::pow(output_radius, degree) * norm /
                    std::pow(
                        input_radius,
                        basis[left].analytic_degree +
                            basis[right].analytic_degree);
                if (score > best_score) {
                    best_score = score;
                    best_degree = degree;
                    best_left = left;
                    best_right = right;
                    best_raw_candidate = std::move(raw_candidate);
                }
            }
        }
        if (!(best_score > entry_tolerance)) {
            break;
        }
        const std::string label =
            "closure-" + std::to_string(extension + 1) + "(" +
            basis[best_left].label + "," +
            basis[best_right].label + ")";
        const int highest_shell =
            LocalSldResponseBasis::highest_active_shell(
                best_raw_candidate);
        std::vector<LocalSldResponseBasisElement> candidates = basis;
        candidates.push_back(LocalSldResponseBasisElement{
            std::move(best_raw_candidate), label, -1,
            highest_shell, best_degree, false});
        basis = LocalSldResponseBasis::graded_orthonormalize(
            dynamics, std::move(candidates));
        ++report.closure_extensions_constructed;
    }
    report.basis_size = basis.size();
    report.basis.reserve(basis.size());
    for (std::size_t index = 0; index < basis.size(); ++index) {
        report.basis.push_back(LocalSldResponseTensorBasisRow{
            static_cast<int>(index),
            basis[index].response_order,
            basis[index].highest_active_shell,
            basis[index].analytic_degree,
            basis[index].label,
            basis[index].scalar_response});
    }
    report.maximum_gram_error =
        LocalSldResponseBasis::maximum_gram_error(basis);
    std::map<int, LocalSldResponseTensorDegreeOffsetRow> offset_envelope;
    report.pairs.reserve(basis.size() * basis.size());
    for (std::size_t left = 0; left < basis.size(); ++left) {
        for (std::size_t right = 0; right < basis.size(); ++right) {
            const SpectralIncrement product =
                dynamics.advection_bilinear_direct_partition(
                    basis.front().state,
                    basis[left].state.velocity,
                    basis[right].state.velocity,
                    TriadPartition::all);
            LocalSldResponseTensorPair pair;
            pair.left_basis_index = static_cast<int>(left);
            pair.right_basis_index = static_cast<int>(right);
            pair.left_order = basis[left].response_order;
            pair.right_order = basis[right].response_order;
            pair.left_analytic_degree = basis[left].analytic_degree;
            pair.right_analytic_degree = basis[right].analytic_degree;
            pair.left_label = basis[left].label;
            pair.right_label = basis[right].label;
            pair.left_highest_shell = basis[left].highest_active_shell;
            pair.right_highest_shell = basis[right].highest_active_shell;
            const SpectralReal pair_input_weight = std::pow(
                input_radius,
                pair.left_analytic_degree +
                    pair.right_analytic_degree);
            const SpectralReal full_norm2 = std::max(
                0.0L,
                LocalSldCyclicBasis::pairing(product, product));
            SpectralReal projected_norm2 = 0.0L;
            SpectralReal weighted_sum = 0.0L;
            int maximum_degree = 0;
            for (const LocalSldResponseBasisElement& element : basis) {
                maximum_degree = std::max(
                    maximum_degree, element.analytic_degree);
            }
            std::vector<SpectralReal> degree_norm2(
                static_cast<std::size_t>(maximum_degree + 1), 0.0L);
            std::map<int, SpectralReal> offset_weighted_norm2;
            SpectralIncrement complement = product;
            for (std::size_t output_index = 0;
                 output_index < basis.size(); ++output_index) {
                const SpectralReal coefficient =
                    LocalSldCyclicBasis::pairing(
                        basis[output_index].state.velocity,
                        product);
                projected_norm2 += coefficient * coefficient;
                weighted_sum += std::pow(
                    output_radius,
                    basis[output_index].analytic_degree) *
                    std::abs(coefficient);
                degree_norm2[static_cast<std::size_t>(
                    basis[output_index].analytic_degree)] +=
                    coefficient * coefficient;
                if (std::abs(coefficient) > entry_tolerance) {
                    const int degree_excess =
                        basis[output_index].analytic_degree -
                        (pair.left_analytic_degree +
                         pair.right_analytic_degree + 1);
                    LocalSldResponseTensorDegreeOffsetRow& offset_row =
                        offset_envelope[degree_excess];
                    offset_row.degree_offset = degree_excess;
                    ++offset_row.retained_entries;
                    offset_row.maximum_absolute_coefficient = std::max(
                        offset_row.maximum_absolute_coefficient,
                        std::abs(coefficient));
                    const SpectralReal weighted_entry = std::pow(
                        output_radius,
                        basis[output_index].analytic_degree) *
                        std::abs(coefficient) / pair_input_weight;
                    offset_row.maximum_weighted_entry = std::max(
                        offset_row.maximum_weighted_entry,
                        weighted_entry);
                    offset_weighted_norm2[degree_excess] +=
                        weighted_entry * weighted_entry;
                    if (degree_excess >
                        pair.maximum_output_degree_excess) {
                        pair.maximum_output_degree_excess = degree_excess;
                        pair.maximum_degree_excess_coefficient =
                            std::abs(coefficient);
                    } else if (degree_excess ==
                               pair.maximum_output_degree_excess &&
                               degree_excess > 0) {
                        pair.maximum_degree_excess_coefficient = std::max(
                            pair.maximum_degree_excess_coefficient,
                            std::abs(coefficient));
                    }
                    pair.entries.push_back(
                        LocalSldResponseTensorEntry{
                            static_cast<int>(output_index),
                            basis[output_index].response_order,
                            basis[output_index].analytic_degree,
                            basis[output_index].label,
                            coefficient});
                    ++report.retained_tensor_entries;
                }
                for (std::size_t mode = 0;
                     mode < complement.size(); ++mode) {
                    for (std::size_t component = 0;
                         component < 3; ++component) {
                        complement[mode][component] -= coefficient *
                            basis[output_index].state
                                .velocity[mode][component];
                    }
                }
            }
            for (const auto& [offset, weighted_norm2] :
                 offset_weighted_norm2) {
                LocalSldResponseTensorDegreeOffsetRow& row =
                    offset_envelope[offset];
                row.maximum_pair_block_contribution = std::max(
                    row.maximum_pair_block_contribution,
                    std::sqrt(std::max(0.0L, weighted_norm2)));
            }
            const SpectralReal complement_norm2 = std::max(
                0.0L,
                LocalSldCyclicBasis::pairing(
                    complement, complement));
            pair.full_norm = std::sqrt(full_norm2);
            pair.projected_norm = std::sqrt(
                std::max(0.0L, projected_norm2));
            pair.complement_norm = std::sqrt(complement_norm2);
            pair.complement_fraction = pair.full_norm > 0.0L
                ? pair.complement_norm / pair.full_norm
                : 0.0L;
            pair.weighted_projected_ratio = weighted_sum /
                std::pow(
                    input_radius,
                    pair.left_analytic_degree +
                        pair.right_analytic_degree);
            SpectralReal weighted_degree_block_sum = 0.0L;
            for (int degree = 0; degree <= maximum_degree; ++degree) {
                weighted_degree_block_sum += std::pow(
                    output_radius, degree) * std::sqrt(std::max(
                        0.0L,
                        degree_norm2[static_cast<std::size_t>(degree)]));
            }
            pair.weighted_projected_degree_block_ratio =
                weighted_degree_block_sum /
                std::pow(
                    input_radius,
                    pair.left_analytic_degree +
                        pair.right_analytic_degree);
            pair.input_weighted_complement_norm =
                pair.complement_norm /
                std::pow(
                    input_radius,
                    pair.left_analytic_degree +
                        pair.right_analytic_degree);
            std::vector<SpectralReal> shell_norm2(
                static_cast<std::size_t>(cutoff + 1), 0.0L);
            SpectralReal h1_norm2 = 0.0L;
            SpectralReal h2_norm2 = 0.0L;
            for (std::size_t mode = 0;
                 mode < complement.size(); ++mode) {
                const WaveVector wave = basis.front().state.waves[mode];
                const int shell = std::max(
                    {std::abs(wave.x), std::abs(wave.y),
                     std::abs(wave.z)});
                const SpectralReal mode_norm2 = std::max(
                    0.0L,
                    std::real(dot_hermitian(
                        complement[mode], complement[mode])));
                shell_norm2[static_cast<std::size_t>(shell)] +=
                    mode_norm2;
                const SpectralReal wave2 = static_cast<SpectralReal>(
                    norm_squared(wave));
                h1_norm2 += wave2 * mode_norm2;
                h2_norm2 += wave2 * wave2 * mode_norm2;
            }
            SpectralReal shell_reconstructed_norm2 = 0.0L;
            SpectralReal shell_weighted_norm = 0.0L;
            for (int shell = 1; shell <= cutoff; ++shell) {
                const SpectralReal norm2 =
                    shell_norm2[static_cast<std::size_t>(shell)];
                shell_reconstructed_norm2 += norm2;
                if (!(norm2 > entry_tolerance * entry_tolerance)) {
                    continue;
                }
                const SpectralReal shell_norm = std::sqrt(norm2);
                const SpectralReal output_weighted_norm =
                    std::pow(output_radius, shell - 1) * shell_norm;
                shell_weighted_norm += output_weighted_norm;
                pair.complement_shells.push_back(
                    LocalSldResponseTensorShell{
                        shell,
                        shell_norm,
                        pair.complement_norm > 0.0L
                            ? shell_norm / pair.complement_norm
                            : 0.0L,
                        output_weighted_norm});
            }
            const SpectralReal input_weight =
                std::pow(
                    input_radius,
                    pair.left_analytic_degree +
                        pair.right_analytic_degree);
            pair.shell_weighted_complement_ratio =
                shell_weighted_norm / input_weight;
            pair.h1_complement_norm = std::sqrt(
                std::max(0.0L, h1_norm2));
            pair.h2_complement_norm = std::sqrt(
                std::max(0.0L, h2_norm2));
            pair.input_weighted_h1_complement_norm =
                pair.h1_complement_norm / input_weight;
            pair.input_weighted_h2_complement_norm =
                pair.h2_complement_norm / input_weight;
            pair.shell_norm_reconstruction_error = std::abs(
                complement_norm2 - shell_reconstructed_norm2) /
                std::max(complement_norm2, 1e-30L);
            const SpectralReal reconstruction_error = std::abs(
                full_norm2 - projected_norm2 - complement_norm2) /
                std::max(full_norm2, 1e-30L);
            report.maximum_norm_reconstruction_error = std::max(
                report.maximum_norm_reconstruction_error,
                reconstruction_error);
            report.maximum_projected_bilinear_constant = std::max(
                report.maximum_projected_bilinear_constant,
                pair.weighted_projected_ratio);
            report.maximum_projected_degree_block_constant = std::max(
                report.maximum_projected_degree_block_constant,
                pair.weighted_projected_degree_block_ratio);
            report.maximum_complement_norm = std::max(
                report.maximum_complement_norm,
                pair.complement_norm);
            report.maximum_complement_fraction = std::max(
                report.maximum_complement_fraction,
                pair.complement_fraction);
            report.maximum_input_weighted_complement_norm = std::max(
                report.maximum_input_weighted_complement_norm,
                pair.input_weighted_complement_norm);
            report.maximum_shell_weighted_complement_ratio = std::max(
                report.maximum_shell_weighted_complement_ratio,
                pair.shell_weighted_complement_ratio);
            report.maximum_input_weighted_h1_complement_norm = std::max(
                report.maximum_input_weighted_h1_complement_norm,
                pair.input_weighted_h1_complement_norm);
            report.maximum_input_weighted_h2_complement_norm = std::max(
                report.maximum_input_weighted_h2_complement_norm,
                pair.input_weighted_h2_complement_norm);
            report.maximum_shell_norm_reconstruction_error = std::max(
                report.maximum_shell_norm_reconstruction_error,
                pair.shell_norm_reconstruction_error);
            if (pair.maximum_output_degree_excess >
                report.maximum_output_degree_excess) {
                report.maximum_output_degree_excess =
                    pair.maximum_output_degree_excess;
                report.maximum_degree_excess_coefficient =
                    pair.maximum_degree_excess_coefficient;
            } else if (pair.maximum_output_degree_excess ==
                           report.maximum_output_degree_excess &&
                       pair.maximum_output_degree_excess > 0) {
                report.maximum_degree_excess_coefficient = std::max(
                    report.maximum_degree_excess_coefficient,
                    pair.maximum_degree_excess_coefficient);
            }
            report.pairs.push_back(std::move(pair));
        }
    }
    report.axis_pair_candidate_constant =
        output_radius / std::sqrt(3.0L);
    report.finite_combined_weighted_bound =
        report.maximum_projected_degree_block_constant +
        report.maximum_shell_weighted_complement_ratio;
    report.projected_constant_over_axis_candidate =
        report.axis_pair_candidate_constant > 0.0L
        ? report.maximum_projected_degree_block_constant /
              report.axis_pair_candidate_constant
        : 0.0L;
    report.axis_pair_candidate_survives =
        report.projected_constant_over_axis_candidate <=
        1.0L + 1e-12L;
    report.graded_support_closed =
        report.maximum_output_degree_excess == 0;
    report.degree_offset_envelope.reserve(offset_envelope.size());
    for (const auto& [offset, row] : offset_envelope) {
        static_cast<void>(offset);
        report.degree_offset_envelope.push_back(row);
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
        } else if (name == "--include-211-transverse") {
            options.include_transverse_two_one_one = true;
        } else if (name == "--include-310-orbits") {
            options.include_three_one_zero_orbits = true;
        } else if (name == "--closure-extensions") {
            options.closure_extensions = std::stoi(next(index, name));
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
        (options.include_three_one_zero_orbits && options.cutoff < 3) ||
        options.closure_extensions < 0 ||
        options.closure_extensions > 16 ||
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
        << "  --include-211-transverse  append transverse (2,1,1) orbit\n"
        << "  --include-310-orbits      append both oriented (3,1,0) orbits\n"
        << "  --closure-extensions N    greedily add 0..16 residual products\n"
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
            options.entry_tolerance,
            options.include_transverse_two_one_one,
            options.include_three_one_zero_orbits,
            options.closure_extensions);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "response tensor cutoff=" << report.cutoff
        << " depth=" << report.depth
        << " radii=" << static_cast<double>(report.input_radius)
        << "->" << static_cast<double>(report.output_radius)
        << " entries=" << report.retained_tensor_entries
        << " basis=" << report.basis_size
        << " closure_extensions="
        << report.closure_extensions_constructed
        << " projected_constant="
        << static_cast<double>(
               report.maximum_projected_bilinear_constant)
        << " projected_block_constant="
        << static_cast<double>(
               report.maximum_projected_degree_block_constant)
        << " axis_candidate="
        << static_cast<double>(report.axis_pair_candidate_constant)
        << " axis_survives="
        << (report.axis_pair_candidate_survives ? "yes" : "no")
        << " complement_norm="
        << static_cast<double>(report.maximum_complement_norm)
        << " complement_fraction="
        << static_cast<double>(report.maximum_complement_fraction)
        << " weighted_complement="
        << static_cast<double>(
               report.maximum_input_weighted_complement_norm)
        << " shell_weighted_complement="
        << static_cast<double>(
               report.maximum_shell_weighted_complement_ratio)
        << " weighted_H1_complement="
        << static_cast<double>(
               report.maximum_input_weighted_h1_complement_norm)
        << " combined_bound="
        << static_cast<double>(report.finite_combined_weighted_bound)
        << " degree_excess="
        << report.maximum_output_degree_excess
        << " Gram_error="
        << static_cast<double>(report.maximum_gram_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
