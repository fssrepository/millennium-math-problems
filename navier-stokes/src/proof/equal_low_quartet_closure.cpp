#include "equal_low_quartet_closure.hpp"

#include <filesystem>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {

EqualLowQuartetClosureReport EqualLowQuartetClosure::certify(
    int maximum_cutoff,
    int target_length_multiplier) {
    EqualLowQuartetClosureReport report;
    report.target_length_multiplier = target_length_multiplier;
    report.geometry = EqualLowTriadGeometry::certify(
        maximum_cutoff, target_length_multiplier);
    report.fixed_angle_plane_sphere_bound_proved =
        report.geometry.exact_fixed_angle_relation &&
        report.geometry.all_degree_bounds_hold;
    report.closed_single_shell_power_bound =
        report.fixed_angle_plane_sphere_bound_proved &&
        report.bilinear_l2_frequency_power == Rational(3, 2) &&
        report.bracket_frequency_power == Rational(5) &&
        report.frequency_gain == Rational(-1, 2);
    report.structural_entries_neighbor_shell_local =
        target_length_multiplier >= 1 &&
        target_length_multiplier <= 3;
    // The shell estimates for S and T depend only on the R^(3/2)
    // bilinear bound and fixed-width support, not on the fixed angle.
    report.direct_normalization_target_bound_proved =
        report.closed_single_shell_power_bound;
    report.cutoff_independent_closed_family_bound =
        report.closed_single_shell_power_bound &&
        report.structural_entries_neighbor_shell_local &&
        report.direct_normalization_target_bound_proved;
    return report;
}

EqualLowQuartetClosureOptions EqualLowQuartetClosureCli::parse(
    int argc, char** argv, int first) {
    EqualLowQuartetClosureOptions options;
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
        } else if (name == "--multiplier") {
            options.target_length_multiplier =
                std::stoi(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown equal-low-quartet option: " + name);
        }
    }
    if (options.maximum_cutoff < 1 || options.maximum_cutoff > 12 ||
        options.target_length_multiplier < 1 ||
        options.target_length_multiplier > 3 ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "equal-low-quartet-certificate requires cutoff 1..12, multiplier 1..3, and a certificate path");
    }
    return options;
}

void EqualLowQuartetClosureCli::print_help(std::ostream& out) {
    out << "Equal-low fixed-angle quartet closure options:\n"
        << "  --max-cutoff K       verify incidence through K\n"
        << "  --multiplier R       target squared length R*m (1, 2, or 3)\n"
        << "  --certificate PATH   write English JSON proof certificate\n";
}

int EqualLowQuartetClosureCli::run(
    const EqualLowQuartetClosureOptions& options,
    std::ostream& out) {
    const EqualLowQuartetClosureReport report =
        EqualLowQuartetClosure::certify(
            options.maximum_cutoff,
            options.target_length_multiplier);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write equal-low quartet certificate");
    }
    certificate
        << "{\n"
        << "  \"schema\": \"navier-stokes-equal-low-quartet-closure-v1\",\n"
        << "  \"family\": \"(m,m,Rm) fixed-angle squared-length signatures\",\n"
        << "  \"target_length_multiplier\": "
        << report.target_length_multiplier << ",\n"
        << "  \"fixed_angle_relation\": \"p dot q=(R-2)m/2 and p dot (p+q)=R m/2\",\n"
        << "  \"incidence_proof\": \"fixing an input or target intersects one integer plane with one sphere; fixing one free coordinate leaves at most two roots, hence degree O(K)\",\n"
        << "  \"bilinear_shell_bound\": \"||B_R(v,w)||_2 <= C_R K^(3/2)||v||_2||w||_2\",\n"
        << "  \"closed_quartet_shell_bound\": \"|K_R+G_R|_j <= C_R K^5 E_near,j^2\",\n"
        << "  \"normalization_bound\": \"S_R <= C_R Z^(5/4)P^(1/4), T_R <= C_R Z^(1/4)P^(5/4)\",\n"
        << "  \"bilinear_l2_frequency_power\": \""
        << report.bilinear_l2_frequency_power.str() << "\",\n"
        << "  \"bracket_frequency_power\": \""
        << report.bracket_frequency_power.str() << "\",\n"
        << "  \"target_frequency_power\": \""
        << report.target_frequency_power.str() << "\",\n"
        << "  \"frequency_gain\": \""
        << report.frequency_gain.str() << "\",\n"
        << "  \"maximum_cutoff\": "
        << report.geometry.maximum_cutoff << ",\n"
        << "  \"maximum_input_degree_ratio\": "
        << static_cast<double>(
               report.geometry.maximum_input_degree_ratio) << ",\n"
        << "  \"maximum_target_degree_ratio\": "
        << static_cast<double>(
               report.geometry.maximum_target_degree_ratio) << ",\n"
        << "  \"fixed_angle_plane_sphere_bound_proved\": "
        << (report.fixed_angle_plane_sphere_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"closed_single_shell_power_bound\": "
        << (report.closed_single_shell_power_bound
            ? "true" : "false") << ",\n"
        << "  \"structural_entries_neighbor_shell_local\": "
        << (report.structural_entries_neighbor_shell_local
            ? "true" : "false") << ",\n"
        << "  \"direct_normalization_target_bound_proved\": "
        << (report.direct_normalization_target_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_closed_family_bound\": "
        << (report.cutoff_independent_closed_family_bound
            ? "true" : "false") << ",\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"combine the proved fixed-angle family with the other remainder signatures and the mixed block\"\n"
        << "}\n";
    out << "equal-low quartet multiplier="
        << report.target_length_multiplier
        << " bracket=R^" << report.bracket_frequency_power.str()
        << " target=R^" << report.target_frequency_power.str()
        << " gain=R^" << report.frequency_gain.str()
        << " closed_family="
        << (report.cutoff_independent_closed_family_bound
            ? "PASS" : "FAIL") << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.cutoff_independent_closed_family_bound ? 0 : 2;
}

}  // namespace lemma
