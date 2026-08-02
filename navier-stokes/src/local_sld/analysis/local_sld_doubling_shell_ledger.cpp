#include "local_sld_doubling_shell_ledger.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"
#include "local_sld_two_scale_state.hpp"

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

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "doubling-shell increment layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

SpectralReal increment_relative_error(
    const SpectralIncrement& reference,
    const SpectralIncrement& reconstructed) {
    SpectralReal difference2 = 0.0L;
    SpectralReal reference2 = 0.0L;
    for (std::size_t mode = 0; mode < reference.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            difference2 += std::norm(
                reference[mode][component] -
                reconstructed[mode][component]);
            reference2 += std::norm(reference[mode][component]);
        }
    }
    return std::sqrt(difference2 /
        std::max(reference2,
                 std::numeric_limits<SpectralReal>::min()));
}

void add_increment(
    SpectralIncrement& target,
    const SpectralIncrement& source) {
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += source[mode][component];
        }
    }
}

void write_certificate(
    const LocalSldDoublingShellReport& report,
    const LocalSldDoublingShellOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write doubling-shell certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-doubling-shell-v2\",\n"
        << "  \"definition\": \"exact ordered dyadic-shell matrix of the complete closed (m,m,2m) K+G bracket\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"evolved_steps\": " << options.evolve_steps << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(options.viscosity) << ",\n"
        << "  \"time_step\": "
        << static_cast<double>(options.time_step) << ",\n"
        << "  \"generated_two_scale\": "
        << (options.two_scale_axis > 0 ? "true" : "false") << ",\n"
        << "  \"two_scale_axis\": " << options.two_scale_axis << ",\n"
        << "  \"high_to_low_energy_ratio\": "
        << static_cast<double>(options.high_to_low_energy_ratio) << ",\n"
        << "  \"two_scale_response_angle\": "
        << static_cast<double>(options.two_scale_response_angle) << ",\n"
        << "  \"enstrophy\": "
        << static_cast<double>(report.enstrophy) << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(report.palinstrophy) << ",\n"
        << "  \"family_stretching\": "
        << static_cast<double>(report.family_stretching) << ",\n"
        << "  \"family_palinstrophy_cross\": "
        << static_cast<double>(report.family_palinstrophy_cross) << ",\n"
        << "  \"absolute_lqc7_ratio\": "
        << static_cast<double>(report.absolute_lqc7_ratio) << ",\n"
        << "  \"signed_local_sld_ratio\": "
        << static_cast<double>(report.signed_local_sld_ratio) << ",\n"
        << "  \"absolute_target_scale_ratio\": "
        << static_cast<double>(report.absolute_target_scale_ratio) << ",\n"
        << "  \"projected_square\": {\n"
        << "    \"expanded_negative_square\": "
        << static_cast<double>(
               report.projected_square.expanded_negative_square) << ",\n"
        << "    \"expanded_enstrophy_remainder\": "
        << static_cast<double>(
               report.projected_square.expanded_enstrophy_remainder)
        << ",\n"
        << "    \"expanded_palinstrophy_remainder\": "
        << static_cast<double>(
               report.projected_square.expanded_palinstrophy_remainder)
        << ",\n"
        << "    \"expanded_total\": "
        << static_cast<double>(report.projected_square.expanded_total)
        << ",\n"
        << "    \"completion_coefficient\": "
        << static_cast<double>(
               report.projected_square.completion_coefficient) << ",\n"
        << "    \"completed_negative_square\": "
        << static_cast<double>(
               report.projected_square.completed_negative_square) << ",\n"
        << "    \"completed_enstrophy_remainder\": "
        << static_cast<double>(
               report.projected_square.completed_enstrophy_remainder)
        << ",\n"
        << "    \"completed_hyperpalinstrophy_remainder\": "
        << static_cast<double>(report.projected_square
                                   .completed_hyperpalinstrophy_remainder)
        << ",\n"
        << "    \"completed_total\": "
        << static_cast<double>(report.projected_square.completed_total)
        << ",\n"
        << "    \"completion_relative_error\": "
        << static_cast<double>(
               report.projected_square.completion_relative_error) << ",\n"
        << "    \"identity_verified\": "
        << (report.projected_square.identity_verified ? "true" : "false")
        << "\n  },\n"
        << "  \"term_totals\": {\n"
        << "    \"outer_square\": "
        << static_cast<double>(report.outer_square_total) << ",\n"
        << "    \"advected_commutator\": "
        << static_cast<double>(report.advected_commutator_total) << ",\n"
        << "    \"enstrophy_normalization\": "
        << static_cast<double>(report.enstrophy_normalization_total)
        << ",\n"
        << "    \"palinstrophy_normalization\": "
        << static_cast<double>(report.palinstrophy_normalization_total)
        << ",\n"
        << "    \"advecting_nested\": "
        << static_cast<double>(report.advecting_nested_total) << "\n"
        << "  },\n"
        << "  \"non_projected_remainder\": "
        << static_cast<double>(report.non_projected_remainder) << ",\n"
        << "  \"non_projected_target_scale_ratio\": "
        << static_cast<double>(report.non_projected_target_scale_ratio)
        << ",\n"
        << "  \"projected_matrix_reconstruction_error\": "
        << static_cast<double>(
               report.projected_matrix_reconstruction_error) << ",\n"
        << "  \"non_projected_reconstruction_error\": "
        << static_cast<double>(
               report.non_projected_reconstruction_error) << ",\n"
        << "  \"full_family_bracket\": "
        << static_cast<double>(report.full_family_bracket) << ",\n"
        << "  \"reconstructed_family_bracket\": "
        << static_cast<double>(report.reconstructed_family_bracket)
        << ",\n"
        << "  \"bracket_reconstruction_error\": "
        << static_cast<double>(report.bracket_reconstruction_error)
        << ",\n"
        << "  \"advection_reconstruction_error\": "
        << static_cast<double>(report.advection_reconstruction_error)
        << ",\n"
        << "  \"diagonal_total\": "
        << static_cast<double>(report.diagonal_total) << ",\n"
        << "  \"off_diagonal_total\": "
        << static_cast<double>(report.off_diagonal_total) << ",\n"
        << "  \"absolute_matrix_sum\": "
        << static_cast<double>(report.absolute_matrix_sum) << ",\n"
        << "  \"signed_cancellation_fraction\": "
        << static_cast<double>(report.signed_cancellation_fraction)
        << ",\n"
        << "  \"maximum_absolute_entry\": "
        << static_cast<double>(report.maximum_absolute_entry) << ",\n"
        << "  \"maximum_absolute_far_structural_entry\": "
        << static_cast<double>(
               report.maximum_absolute_far_structural_entry) << ",\n"
        << "  \"shells\": [\n";
    for (std::size_t index = 0; index < report.shells.size(); ++index) {
        const LocalSldDoublingShellRow& row = report.shells[index];
        output << "    {\"shell\": " << row.shell
            << ", \"minimum_input_squared\": "
            << row.minimum_input_squared
            << ", \"maximum_input_squared_exclusive\": "
            << row.maximum_input_squared_exclusive
            << ", \"advection_norm2\": "
            << static_cast<double>(row.advection_norm2)
            << ", \"stretching\": "
            << static_cast<double>(row.stretching)
            << ", \"palinstrophy_cross\": "
            << static_cast<double>(row.palinstrophy_cross) << '}'
            << (index + 1 == report.shells.size() ? "\n" : ",\n");
    }
    output << "  ],\n  \"matrix\": [\n";
    for (std::size_t index = 0; index < report.matrix.size(); ++index) {
        const LocalSldDoublingShellMatrixEntry& entry =
            report.matrix[index];
        output << "    {\"left_shell\": " << entry.left_shell
            << ", \"right_shell\": " << entry.right_shell
            << ", \"outer_square\": "
            << static_cast<double>(entry.outer_square)
            << ", \"advected_commutator\": "
            << static_cast<double>(entry.advected_commutator)
            << ", \"enstrophy_normalization\": "
            << static_cast<double>(entry.enstrophy_normalization)
            << ", \"palinstrophy_normalization\": "
            << static_cast<double>(entry.palinstrophy_normalization)
            << ", \"advecting_nested\": "
            << static_cast<double>(entry.advecting_nested)
            << ", \"total\": " << static_cast<double>(entry.total)
            << '}'
            << (index + 1 == report.matrix.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"exact_reconstruction\": "
        << (report.exact_reconstruction ? "true" : "false") << ",\n"
        << "  \"structural_entries_are_neighbor_shell_local\": "
        << (report.structural_entries_are_neighbor_shell_local
            ? "true" : "false") << ",\n"
        << "  \"signed_cross_shell_bound_proved\": false,\n"
        << "  \"finite_matrix_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldDoublingShellReport LocalSldDoublingShellLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state) {
    LocalSldDoublingShellReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    const TriadSelection full_selection =
        TriadSelection::local_equal_low_doubling();
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics, full_selection).evaluate(state);
    report.enstrophy = full.enstrophy;
    report.palinstrophy = full.palinstrophy;
    report.family_stretching = full.signed_stretching;
    report.family_palinstrophy_cross = full.palinstrophy_cross;
    report.absolute_lqc7_ratio = full.constant_ratio;
    report.signed_local_sld_ratio = full.signed_local_sld_ratio;
    report.full_family_bracket = full.signed_two_entry_bracket;
    const SpectralReal target_scale =
        std::pow(report.enstrophy, 1.25L) *
        std::pow(report.palinstrophy, 0.75L);
    if (target_scale > 0.0L) {
        report.absolute_target_scale_ratio =
            std::abs(report.full_family_bracket) / target_scale;
    }
    report.projected_square = LocalSldProjectedSquare::evaluate(
        dynamics, state, full_selection);
    report.non_projected_remainder = report.full_family_bracket -
        report.projected_square.expanded_total;
    if (target_scale > 0.0L) {
        report.non_projected_target_scale_ratio =
            std::abs(report.non_projected_remainder) / target_scale;
    }
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const SpectralIncrement full_advection =
        dynamics.advection_direct_partition(state, full_selection);
    SpectralIncrement reconstructed_advection(
        state.velocity.size(), ComplexVector{});

    struct ShellWork {
        TriadSelection selection;
        SpectralIncrement advection;
        SpectralIncrement weighted_advection;
        SpectralIncrement transported_au;
        SpectralReal stretching = 0.0L;
        SpectralReal palinstrophy_cross = 0.0L;
    };
    std::vector<ShellWork> work;
    SpectralInteger lower = 1;
    const SpectralInteger maximum_squared =
        3 * static_cast<SpectralInteger>(report.cutoff) * report.cutoff;
    for (int shell = 0; lower <= maximum_squared; ++shell) {
        const SpectralInteger upper = lower >
                std::numeric_limits<SpectralInteger>::max() / 4
            ? std::numeric_limits<SpectralInteger>::max()
            : 4 * lower;
        const TriadSelection selection =
            TriadSelection::local_equal_low_doubling_shell(
                lower, upper);
        ShellWork shell_work{
            selection,
            dynamics.advection_direct_partition(state, selection),
            {}, {}, 0.0L, 0.0L};
        shell_work.weighted_advection = laplacian_weight(
            state, shell_work.advection);
        shell_work.transported_au =
            dynamics.advection_bilinear_direct_partition(
                state, state.velocity, au, selection);
        shell_work.stretching = pairing(au, shell_work.advection);
        shell_work.palinstrophy_cross = pairing(
            shell_work.weighted_advection, au);
        add_increment(reconstructed_advection, shell_work.advection);
        report.shells.push_back(LocalSldDoublingShellRow{
            shell, lower, upper,
            pairing(shell_work.advection, shell_work.advection),
            shell_work.stretching,
            shell_work.palinstrophy_cross});
        work.push_back(std::move(shell_work));
        if (upper == std::numeric_limits<SpectralInteger>::max()) {
            break;
        }
        lower = upper;
    }
    report.advection_reconstruction_error = increment_relative_error(
        full_advection, reconstructed_advection);

    for (std::size_t left = 0; left < work.size(); ++left) {
        for (std::size_t right = 0; right < work.size(); ++right) {
            LocalSldDoublingShellMatrixEntry entry;
            entry.left_shell = static_cast<int>(left);
            entry.right_shell = static_cast<int>(right);
            entry.outer_square = -pairing(
                work[left].advection,
                work[right].weighted_advection);
            entry.advected_commutator = pairing(
                work[left].advection,
                work[right].transported_au);
            entry.enstrophy_normalization =
                work[left].stretching * work[right].stretching /
                (2.0L * report.enstrophy);
            entry.palinstrophy_normalization =
                3.0L * work[left].stretching *
                work[right].palinstrophy_cross /
                (2.0L * report.palinstrophy);
            const SpectralIncrement nested =
                dynamics.advection_bilinear_direct_partition(
                    state,
                    work[left].advection,
                    state.velocity,
                    work[right].selection);
            entry.advecting_nested = -pairing(au, nested);
            entry.total = entry.outer_square +
                entry.advected_commutator +
                entry.enstrophy_normalization +
                entry.palinstrophy_normalization +
                entry.advecting_nested;
            report.outer_square_total += entry.outer_square;
            report.advected_commutator_total +=
                entry.advected_commutator;
            report.enstrophy_normalization_total +=
                entry.enstrophy_normalization;
            report.palinstrophy_normalization_total +=
                entry.palinstrophy_normalization;
            report.advecting_nested_total += entry.advecting_nested;
            report.reconstructed_family_bracket += entry.total;
            report.absolute_matrix_sum += std::abs(entry.total);
            report.maximum_absolute_entry = std::max(
                report.maximum_absolute_entry, std::abs(entry.total));
            if (left == right) {
                report.diagonal_total += entry.total;
            } else {
                report.off_diagonal_total += entry.total;
                const std::size_t gap = left > right
                    ? left - right : right - left;
                if (gap > 1) {
                    report.maximum_absolute_far_structural_entry =
                        std::max({report
                                      .maximum_absolute_far_structural_entry,
                                  std::abs(entry.outer_square),
                                  std::abs(entry.advected_commutator),
                                  std::abs(entry.advecting_nested)});
                }
            }
            report.matrix.push_back(entry);
        }
    }
    report.bracket_reconstruction_error = std::abs(
        report.full_family_bracket -
        report.reconstructed_family_bracket) /
        std::max(std::abs(report.full_family_bracket), 1e-30L);
    const SpectralReal projected_from_matrix =
        report.outer_square_total +
        report.enstrophy_normalization_total +
        report.palinstrophy_normalization_total;
    report.projected_matrix_reconstruction_error = std::abs(
        projected_from_matrix - report.projected_square.expanded_total) /
        std::max({std::abs(projected_from_matrix),
                  std::abs(report.projected_square.expanded_total),
                  1e-30L});
    const SpectralReal non_projected_from_matrix =
        report.advected_commutator_total +
        report.advecting_nested_total;
    report.non_projected_reconstruction_error = std::abs(
        non_projected_from_matrix - report.non_projected_remainder) /
        std::max({std::abs(non_projected_from_matrix),
                  std::abs(report.non_projected_remainder), 1e-30L});
    if (report.absolute_matrix_sum > 0.0L) {
        report.signed_cancellation_fraction = 1.0L -
            std::abs(report.reconstructed_family_bracket) /
                report.absolute_matrix_sum;
    }
    report.exact_reconstruction =
        report.advection_reconstruction_error < 1e-14L &&
        report.bracket_reconstruction_error < 1e-14L &&
        report.projected_matrix_reconstruction_error < 1e-14L &&
        report.non_projected_reconstruction_error < 1e-14L;
    report.structural_entries_are_neighbor_shell_local =
        report.maximum_absolute_far_structural_entry < 1e-14L;
    return report;
}

LocalSldDoublingShellOptions LocalSldDoublingShellCli::parse(
    int argc, char** argv, int first) {
    LocalSldDoublingShellOptions options;
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
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--evolve-steps") {
            options.evolve_steps = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else if (name == "--backend") {
            options.backend = next(index, name);
        } else if (name == "--two-scale-axis") {
            options.two_scale_axis = std::stoi(next(index, name));
        } else if (name == "--high-energy-ratio") {
            options.high_to_low_energy_ratio = std::stold(
                next(index, name));
        } else if (name == "--two-scale-angle") {
            options.two_scale_response_angle = std::stold(
                next(index, name));
        } else if (name == "--state-output") {
            options.state_output_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-sld-doubling-shell option: " + name);
        }
    }
    const bool generated = options.two_scale_axis > 0;
    if (generated && !(options.high_to_low_energy_ratio > 0.0L)) {
        options.high_to_low_energy_ratio = std::pow(
            static_cast<SpectralReal>(options.two_scale_axis),
            -11.0L / 4.0L);
    }
    if ((!generated && options.state_path.empty()) ||
        (generated && (options.two_scale_axis < 2 ||
                       options.two_scale_axis > 12 ||
                       !(options.high_to_low_energy_ratio > 0.0L) ||
                       !std::isfinite(
                           options.high_to_low_energy_ratio) ||
                       !std::isfinite(
                           options.two_scale_response_angle))) ||
        options.certificate_path.empty() ||
        options.threads < 1 || options.threads > 256 ||
        options.evolve_steps < 0 || !(options.viscosity > 0.0L) ||
        !(options.time_step > 0.0L) ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft")) {
        throw std::invalid_argument(
            "local-sld-doubling-shells requires state, certificate, and valid numerical options");
    }
    return options;
}

void LocalSldDoublingShellCli::print_help(std::ostream& out) {
    out << "Local SLD doubling-family shell-matrix options:\n"
        << "  --state PATH         input Fourier TSV\n"
        << "  --threads N          direct/evolution workers\n"
        << "  --evolve-steps N     optional RK4 evolution steps\n"
        << "  --nu X               optional evolution viscosity\n"
        << "  --dt X               optional evolution time step\n"
        << "  --backend NAME       direct, fft, or auto\n"
        << "  --two-scale-axis L   generate cyclic axes at scales 1 and L\n"
        << "  --high-energy-ratio X  generated high/low energy ratio; default L^(-11/4)\n"
        << "  --two-scale-angle X  axis/response mixture angle; default 0.24\n"
        << "  --state-output PATH  optionally save generated state TSV\n"
        << "  --certificate PATH   write English JSON shell matrix\n";
}

int LocalSldDoublingShellCli::run(
    const LocalSldDoublingShellOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure(options.backend, options.threads);
    const SpectralDynamics dynamics(galerkin);
    SpectralState state = options.two_scale_axis > 0
        ? LocalSldTwoScaleState::cyclic_response_mixture(
              dynamics,
              options.two_scale_axis,
              options.high_to_low_energy_ratio,
              options.two_scale_response_angle)
        : SpectralStateReader::read_tsv(options.state_path);
    if (!options.state_output_path.empty()) {
        SpectralStateWriter::write_tsv(
            options.state_output_path, state,
            "two-scale cyclic axis/response state; generated by local-sld-doubling-shells");
    }
    for (int step = 0; step < options.evolve_steps; ++step) {
        dynamics.rk4_step(
            state, options.viscosity, options.time_step);
    }
    const LocalSldDoublingShellReport report =
        LocalSldDoublingShellLedger::analyze(dynamics, state);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "doubling shell matrix cutoff=" << report.cutoff
        << " shells=" << report.shells.size()
        << " bracket=" << static_cast<double>(report.full_family_bracket)
        << " diagonal=" << static_cast<double>(report.diagonal_total)
        << " off_diagonal="
        << static_cast<double>(report.off_diagonal_total)
        << " cancellation="
        << static_cast<double>(report.signed_cancellation_fraction)
        << " target_ratio="
        << static_cast<double>(report.absolute_target_scale_ratio)
        << " signed_sld="
        << static_cast<double>(report.signed_local_sld_ratio)
        << " projected="
        << static_cast<double>(report.projected_square.expanded_total)
        << " remainder="
        << static_cast<double>(report.non_projected_remainder)
        << " reconstruction="
        << static_cast<double>(report.bracket_reconstruction_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_reconstruction ? 0 : 2;
}

}  // namespace lemma
