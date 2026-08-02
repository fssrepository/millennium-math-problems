#include "projective_quartet_closure.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

std::array<SpectralInteger, 3> parse_signature(
    const std::string& text) {
    std::array<SpectralInteger, 3> result{};
    std::stringstream stream(text);
    std::string token;
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (!std::getline(stream, token, ',')) {
            throw std::invalid_argument(
                "projective signature must have three comma-separated integers");
        }
        result[index] = static_cast<SpectralInteger>(std::stoll(token));
    }
    if (std::getline(stream, token, ',')) {
        throw std::invalid_argument(
            "projective signature must have exactly three entries");
    }
    std::sort(result.begin(), result.end());
    return result;
}

}  // namespace

ProjectiveQuartetClosureReport ProjectiveQuartetClosure::certify(
    int maximum_cutoff,
    std::array<SpectralInteger, 3> primitive_squared_lengths) {
    ProjectiveQuartetClosureReport report;
    report.geometry = ProjectiveTriadGeometry::certify(
        maximum_cutoff, primitive_squared_lengths);
    report.fixed_projective_ray_incidence_proved =
        report.geometry.primitive_signature &&
        report.geometry.triangle_feasible &&
        report.geometry.fixed_plane_sphere_geometry &&
        report.geometry.all_degree_bounds_hold;
    report.fixed_projective_ray_shell_bound_proved =
        report.fixed_projective_ray_incidence_proved &&
        report.incidence_degree_power == Rational(1) &&
        report.bilinear_l2_frequency_power == Rational(3, 2) &&
        report.bracket_frequency_power == Rational(5) &&
        report.frequency_gain == Rational(-1, 2);
    report.fixed_projective_ray_normalization_bound_proved =
        report.fixed_projective_ray_shell_bound_proved;
    report.cutoff_independent_fixed_projective_ray_bound =
        report.fixed_projective_ray_shell_bound_proved &&
        report.fixed_projective_ray_normalization_bound_proved;
    return report;
}

ProjectiveQuartetClosureOptions ProjectiveQuartetClosureCli::parse(
    int argc, char** argv, int first) {
    ProjectiveQuartetClosureOptions options;
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
        } else if (name == "--signature") {
            options.primitive_squared_lengths = parse_signature(
                next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown projective-quartet option: " + name);
        }
    }
    if (options.maximum_cutoff < 1 || options.maximum_cutoff > 12 ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "projective-quartet-certificate requires cutoff 1..12 and a certificate path");
    }
    return options;
}

void ProjectiveQuartetClosureCli::print_help(std::ostream& out) {
    out << "Fixed projective-signature quartet closure options:\n"
        << "  --max-cutoff K       verify incidence through K\n"
        << "  --signature A,B,C    primitive sorted squared-length ratios\n"
        << "  --certificate PATH   write English JSON proof certificate\n";
}

int ProjectiveQuartetClosureCli::run(
    const ProjectiveQuartetClosureOptions& options,
    std::ostream& out) {
    const ProjectiveQuartetClosureReport report =
        ProjectiveQuartetClosure::certify(
            options.maximum_cutoff,
            options.primitive_squared_lengths);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective quartet certificate");
    }
    const auto& signature = report.geometry.primitive_squared_lengths;
    certificate
        << "{\n"
        << "  \"schema\": \"navier-stokes-projective-quartet-closure-v1\",\n"
        << "  \"primitive_squared_lengths\": ["
        << signature[0] << ", " << signature[1] << ", "
        << signature[2] << "],\n"
        << "  \"family\": \"(a m,b m,c m) for every positive integer scale m\",\n"
        << "  \"incidence_proof\": \"for every ordered role assignment, fixing an input or target leaves one plane-sphere intersection; fixing one coordinate leaves at most two roots, so the degree is O_shape(K)\",\n"
        << "  \"bilinear_shell_bound\": \"||B_shape(v,w)||_2 <= C_shape K^(3/2)||v||_2||w||_2\",\n"
        << "  \"closed_quartet_shell_bound\": \"|K_shape+G_shape|_j <= C_shape K^5 E_near,j^2\",\n"
        << "  \"normalization_bound\": \"S_shape <= C_shape Z^(5/4)P^(1/4), T_shape <= C_shape Z^(1/4)P^(5/4)\",\n"
        << "  \"maximum_cutoff\": "
        << report.geometry.maximum_cutoff << ",\n"
        << "  \"ordered_role_count\": "
        << report.geometry.ordered_role_count << ",\n"
        << "  \"maximum_input_degree_ratio\": "
        << static_cast<double>(
               report.geometry.maximum_input_degree_ratio) << ",\n"
        << "  \"maximum_target_degree_ratio\": "
        << static_cast<double>(
               report.geometry.maximum_target_degree_ratio) << ",\n"
        << "  \"incidence_degree_power\": \""
        << report.incidence_degree_power.str() << "\",\n"
        << "  \"bilinear_l2_frequency_power\": \""
        << report.bilinear_l2_frequency_power.str() << "\",\n"
        << "  \"bracket_frequency_power\": \""
        << report.bracket_frequency_power.str() << "\",\n"
        << "  \"target_frequency_power\": \""
        << report.target_frequency_power.str() << "\",\n"
        << "  \"frequency_gain\": \""
        << report.frequency_gain.str() << "\",\n"
        << "  \"fixed_projective_ray_incidence_proved\": "
        << (report.fixed_projective_ray_incidence_proved
            ? "true" : "false") << ",\n"
        << "  \"fixed_projective_ray_shell_bound_proved\": "
        << (report.fixed_projective_ray_shell_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"fixed_projective_ray_normalization_bound_proved\": "
        << (report.fixed_projective_ray_normalization_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_fixed_projective_ray_bound\": "
        << (report.cutoff_independent_fixed_projective_ray_bound
            ? "true" : "false") << ",\n"
        << "  \"uniform_sum_over_projective_shapes_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"derive shape-uniform constants and a summable estimate over all primitive projective signatures\"\n"
        << "}\n";
    out << "projective quartet signature="
        << signature[0] << ',' << signature[1] << ',' << signature[2]
        << " bracket=R^" << report.bracket_frequency_power.str()
        << " target=R^" << report.target_frequency_power.str()
        << " fixed_ray="
        << (report.cutoff_independent_fixed_projective_ray_bound
            ? "PASS" : "FAIL") << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.cutoff_independent_fixed_projective_ray_bound ? 0 : 2;
}

}  // namespace lemma
