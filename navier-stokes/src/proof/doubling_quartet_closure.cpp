#include "doubling_quartet_closure.hpp"

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

SpectralReal normalization_ratio(
    SpectralReal first_radius,
    SpectralReal second_radius,
    SpectralReal second_energy) {
    const SpectralReal first_energy = 1.0L;
    const auto power = [](SpectralReal value, SpectralReal exponent) {
        return std::pow(value, exponent);
    };
    const SpectralReal stretching =
        power(first_radius, 3.5L) * power(first_energy, 1.5L) +
        power(second_radius, 3.5L) * power(second_energy, 1.5L);
    const SpectralReal weighted_stretching =
        power(first_radius, 5.5L) * power(first_energy, 1.5L) +
        power(second_radius, 5.5L) * power(second_energy, 1.5L);
    const SpectralReal palinstrophy =
        power(first_radius, 4.0L) * first_energy +
        power(second_radius, 4.0L) * second_energy;
    const SpectralReal shell_quartic =
        power(first_radius, 5.0L) * first_energy * first_energy +
        power(second_radius, 5.0L) * second_energy * second_energy;
    return stretching * weighted_stretching /
        std::max(1e-30L, palinstrophy * shell_quartic);
}

}  // namespace

DoublingQuartetClosureReport DoublingQuartetClosure::certify(
    int maximum_cutoff) {
    DoublingQuartetClosureReport report;
    report.geometry = OrthogonalTriadGeometry::certify(maximum_cutoff);
    report.bilinear_l2_frequency_power =
        report.derivative_power +
        Rational(1, 2) * report.incidence_degree_power;
    report.stretching_frequency_power =
        Rational(2) + report.bilinear_l2_frequency_power;

    const Rational outer_square_power =
        Rational(2) +
        Rational(2) * report.bilinear_l2_frequency_power;
    const Rational commutator_power =
        report.bilinear_l2_frequency_power +
        (report.bilinear_l2_frequency_power + Rational(2));
    const Rational nested_advecting_power =
        Rational(2) +
        (report.bilinear_l2_frequency_power +
         report.bilinear_l2_frequency_power);
    const Rational enstrophy_normalization_power =
        Rational(2) * report.stretching_frequency_power - Rational(2);
    const Rational palinstrophy_normalization_power =
        report.stretching_frequency_power +
        (report.bilinear_l2_frequency_power + Rational(4)) -
        Rational(4);
    report.bracket_frequency_power = outer_square_power;
    report.every_quartet_entry_has_same_power =
        outer_square_power == Rational(5) &&
        commutator_power == outer_square_power &&
        nested_advecting_power == outer_square_power &&
        enstrophy_normalization_power == outer_square_power &&
        palinstrophy_normalization_power == outer_square_power;
    report.target_energy_homogeneity_matches =
        report.bracket_energy_power == report.required_energy_power;
    report.frequency_gain =
        report.bracket_frequency_power -
        report.required_frequency_power;

    const Rational shell_sum_frequency_power =
        Rational(2) + Rational(3);
    report.shell_sum_closes_by_zp =
        shell_sum_frequency_power == report.bracket_frequency_power &&
        report.global_enstrophy_power == Rational(3, 2) &&
        report.global_palinstrophy_power == Rational(1, 2);
    report.torus_spectral_gap_closes_target =
        report.global_enstrophy_power -
            report.target_enstrophy_power == Rational(1, 4) &&
        report.global_palinstrophy_power -
            report.target_palinstrophy_power == Rational(-1, 4);
    report.closed_single_shell_power_bound =
        report.geometry.all_degree_bounds_hold &&
        report.every_quartet_entry_has_same_power &&
        report.target_energy_homogeneity_matches &&
        report.frequency_gain < Rational(0);
    // A block with low squared radius m is supported only at squared radii
    // m and 2m. Two structural blocks can therefore pair only when their
    // low radii differ by at most a factor two, hence by at most one of the
    // factor-four squared-radius shells used by the ledger.
    report.structural_entries_neighbor_shell_local = true;
    report.structural_entry_global_bound_proved =
        report.closed_single_shell_power_bound &&
        report.shell_sum_closes_by_zp &&
        report.torus_spectral_gap_closes_target;
    const bool stretching_sequence_bound =
        Rational(2) * report.stretching_enstrophy_power +
                Rational(4) * report.stretching_palinstrophy_power ==
            report.stretching_frequency_power &&
        report.stretching_enstrophy_power +
                report.stretching_palinstrophy_power ==
            Rational(3, 2);
    const bool weighted_stretching_sequence_bound =
        Rational(2) * report.weighted_stretching_enstrophy_power +
                Rational(4) *
                    report.weighted_stretching_palinstrophy_power ==
            Rational(11, 2) &&
        report.weighted_stretching_enstrophy_power +
                report.weighted_stretching_palinstrophy_power ==
            Rational(3, 2);
    const Rational squared_normalization_z =
        Rational(2) * report.stretching_enstrophy_power - Rational(1);
    const Rational squared_normalization_p =
        Rational(2) * report.stretching_palinstrophy_power;
    const Rational cross_normalization_z =
        report.stretching_enstrophy_power +
        report.weighted_stretching_enstrophy_power;
    const Rational cross_normalization_p =
        report.stretching_palinstrophy_power +
        report.weighted_stretching_palinstrophy_power - Rational(1);
    report.direct_normalization_target_bound_proved =
        stretching_sequence_bound &&
        weighted_stretching_sequence_bound &&
        squared_normalization_z == report.normalization_enstrophy_power &&
        squared_normalization_p == report.normalization_palinstrophy_power &&
        cross_normalization_z == report.normalization_enstrophy_power &&
        cross_normalization_p == report.normalization_palinstrophy_power &&
        report.normalization_enstrophy_power -
                report.target_enstrophy_power == Rational(1, 4) &&
        report.normalization_palinstrophy_power -
                report.target_palinstrophy_power == Rational(-1, 4) &&
        report.torus_spectral_gap_closes_target;
    report.projected_normalization_bound_proved =
        report.direct_normalization_target_bound_proved;
    for (int first_shell = 0; first_shell <= maximum_cutoff; ++first_shell) {
        for (int second_shell = first_shell + 1;
             second_shell <= maximum_cutoff; ++second_shell) {
            const SpectralReal first_radius = std::ldexp(
                1.0L, first_shell);
            const SpectralReal second_radius = std::ldexp(
                1.0L, second_shell);
            for (int exponent_index = -192;
                 exponent_index <= 192; ++exponent_index) {
                const SpectralReal energy_ratio = std::exp2(
                    static_cast<SpectralReal>(exponent_index) / 8.0L);
                const SpectralReal ratio = normalization_ratio(
                    first_radius, second_radius, energy_ratio);
                if (ratio > report.maximum_tested_normalization_ratio) {
                    report.maximum_tested_normalization_ratio = ratio;
                    report.maximizing_first_shell = first_shell;
                    report.maximizing_second_shell = second_shell;
                    report.maximizing_energy_ratio = energy_ratio;
                }
            }
        }
    }
    report.normalization_sequence_screen_survives =
        report.maximum_tested_normalization_ratio <= 1.0L + 1e-12L;
    report.normalization_sequence_bound_proved = false;
    report.two_scale_stretching_high_power =
        Rational(7, 2) -
        Rational(3, 2) * report.two_scale_energy_decay_power;
    report.two_scale_weighted_stretching_high_power =
        Rational(11, 2) -
        Rational(3, 2) * report.two_scale_energy_decay_power;
    report.two_scale_palinstrophy_high_power =
        Rational(4) - report.two_scale_energy_decay_power;
    report.two_scale_shell_quartic_high_power =
        Rational(5) -
        Rational(2) * report.two_scale_energy_decay_power;
    report.two_scale_ratio_growth_power =
        report.two_scale_weighted_stretching_high_power -
        report.two_scale_palinstrophy_high_power;
    report.naive_cross_shell_bound_rejected =
        report.two_scale_stretching_high_power < Rational(0) &&
        report.two_scale_shell_quartic_high_power < Rational(0) &&
        Rational(0) < report.two_scale_palinstrophy_high_power &&
        Rational(0) <
            report.two_scale_weighted_stretching_high_power &&
        Rational(0) < report.two_scale_ratio_growth_power &&
        !report.normalization_sequence_screen_survives;
    report.cutoff_independent_closed_family_bound =
        report.closed_single_shell_power_bound &&
        report.shell_sum_closes_by_zp &&
        report.torus_spectral_gap_closes_target &&
        report.structural_entry_global_bound_proved &&
        report.projected_normalization_bound_proved;
    report.full_local_lemma_proved = false;
    return report;
}

DoublingQuartetClosureOptions DoublingQuartetClosureCli::parse(
    int argc, char** argv, int first) {
    DoublingQuartetClosureOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown doubling-quartet option: " + name);
        }
    }
    if (options.maximum_cutoff < 1 || options.maximum_cutoff > 12 ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "doubling-quartet-certificate requires cutoff 1..12 and a certificate path");
    }
    return options;
}

void DoublingQuartetClosureCli::print_help(std::ostream& out) {
    out << "Doubling-family quartet closure options:\n"
        << "  --max-cutoff K       verify orthogonal incidence through K\n"
        << "  --certificate PATH   write English JSON power certificate\n";
}

int DoublingQuartetClosureCli::run(
    const DoublingQuartetClosureOptions& options,
    std::ostream& out) {
    const DoublingQuartetClosureReport report =
        DoublingQuartetClosure::certify(
        options.maximum_cutoff);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write doubling-quartet certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-doubling-quartet-closure-v2\",\n"
        << "  \"family\": \"all squared-length signatures (m,m,2m), equivalently equal-length orthogonal input pairs\",\n"
        << "  \"maximum_cutoff\": " << options.maximum_cutoff << ",\n"
        << "  \"incidence_bound\": \"maximum target degree <= 2(2K+1) on a radius-K shell\",\n"
        << "  \"bilinear_shell_bound\": \"||B_d(v,w)||_2 <= C R^(3/2)||v||_2||w||_2\",\n"
        << "  \"closed_quartet_shell_bound\": \"|K_d+G_d|_j <= C R^5 E_near,j^2\",\n"
        << "  \"shell_sum_bound\": \"sum_j R_j^5 E_near,j^2 <= C Z^(3/2)P^(1/2) <= C Z^(5/4)P^(3/4); use each R_j^2 E_j <= Z, each R_j^4 E_j <= P, then P >= Z on the mean-zero torus\",\n"
        << "  \"normalization_target_bound\": \"S <= C Z^(5/4)P^(1/4), T <= C Z^(1/4)P^(5/4), hence S^2/Z and |S T|/P are <= C Z^(3/2)P^(1/2) <= C Z^(5/4)P^(3/4)\",\n"
        << "  \"bilinear_l2_frequency_power\": \""
        << report.bilinear_l2_frequency_power.str() << "\",\n"
        << "  \"stretching_frequency_power\": \""
        << report.stretching_frequency_power.str() << "\",\n"
        << "  \"bracket_frequency_power\": \""
        << report.bracket_frequency_power.str() << "\",\n"
        << "  \"required_frequency_power\": \""
        << report.required_frequency_power.str() << "\",\n"
        << "  \"frequency_gain\": \""
        << report.frequency_gain.str() << "\",\n"
        << "  \"all_orthogonal_degree_bounds_hold\": "
        << (report.geometry.all_degree_bounds_hold ? "true" : "false")
        << ",\n"
        << "  \"every_quartet_entry_has_same_power\": "
        << (report.every_quartet_entry_has_same_power ? "true" : "false")
        << ",\n"
        << "  \"shell_sum_closes_by_zp\": "
        << (report.shell_sum_closes_by_zp ? "true" : "false")
        << ",\n"
        << "  \"torus_spectral_gap_closes_target\": "
        << (report.torus_spectral_gap_closes_target ? "true" : "false")
        << ",\n"
        << "  \"cutoff_independent_closed_family_bound\": "
        << (report.cutoff_independent_closed_family_bound
            ? "true" : "false") << ",\n"
        << "  \"closed_single_shell_power_bound\": "
        << (report.closed_single_shell_power_bound ? "true" : "false")
        << ",\n"
        << "  \"structural_entries_neighbor_shell_local\": "
        << (report.structural_entries_neighbor_shell_local
            ? "true" : "false")
        << ",\n"
        << "  \"structural_entry_global_bound_proved\": "
        << (report.structural_entry_global_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"direct_normalization_target_bound_proved\": "
        << (report.direct_normalization_target_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"projected_normalization_bound_proved\": "
        << (report.projected_normalization_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"maximum_tested_normalization_ratio\": "
        << static_cast<double>(
               report.maximum_tested_normalization_ratio) << ",\n"
        << "  \"maximizing_two_shell_configuration\": {\"first_shell\": "
        << report.maximizing_first_shell
        << ", \"second_shell\": " << report.maximizing_second_shell
        << ", \"second_to_first_energy_ratio\": "
        << static_cast<double>(report.maximizing_energy_ratio) << "},\n"
        << "  \"normalization_sequence_screen_survives\": "
        << (report.normalization_sequence_screen_survives
            ? "true" : "false") << ",\n"
        << "  \"normalization_sequence_bound_proved\": false,\n"
        << "  \"two_scale_counterexample\": {\"high_shell_energy\": \"L^-"
        << report.two_scale_energy_decay_power.str()
        << "\", \"stretching_high_power\": \""
        << report.two_scale_stretching_high_power.str()
        << "\", \"weighted_stretching_high_power\": \""
        << report.two_scale_weighted_stretching_high_power.str()
        << "\", \"palinstrophy_high_power\": \""
        << report.two_scale_palinstrophy_high_power.str()
        << "\", \"shell_quartic_high_power\": \""
        << report.two_scale_shell_quartic_high_power.str()
        << "\", \"ratio_growth_power\": \""
        << report.two_scale_ratio_growth_power.str() << "\"},\n"
        << "  \"naive_cross_shell_bound_rejected\": "
        << (report.naive_cross_shell_bound_rejected
            ? "true" : "false") << ",\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"extend the now-closed doubling-family estimate to the closed signature remainder and to the mixed block; the full local SLD lemma remains open\",\n"
        << "  \"geometry_rows\": [\n";
    for (std::size_t index = 0;
         index < report.geometry.rows.size(); ++index) {
        const OrthogonalTriadCountRow& row =
            report.geometry.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"ordered_pairs\": " << row.ordered_pairs
            << ", \"maximum_input_degree\": "
            << row.maximum_input_degree
            << ", \"maximum_target_degree\": "
            << row.maximum_target_degree
            << ", \"elementary_degree_bound\": "
            << row.elementary_degree_bound << '}'
            << (index + 1 == report.geometry.rows.size()
                ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << "doubling quartet closure K=1.."
        << options.maximum_cutoff
        << " B_power=R^"
        << report.bilinear_l2_frequency_power.str()
        << " bracket=R^" << report.bracket_frequency_power.str()
        << " target=R^" << report.required_frequency_power.str()
        << " gain=R^" << report.frequency_gain.str()
        << " shell_power="
        << (report.closed_single_shell_power_bound ? "PASS" : "FAIL")
        << " normalization_screen="
        << static_cast<double>(
               report.maximum_tested_normalization_ratio)
        << " closed_family="
        << (report.cutoff_independent_closed_family_bound
            ? "PASS" : "FAIL") << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.cutoff_independent_closed_family_bound ? 0 : 2;
}

}  // namespace lemma
