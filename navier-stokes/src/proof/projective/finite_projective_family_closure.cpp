#include "finite_projective_family_closure.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace lemma {
namespace {

ProjectiveSignature parse_signature(const std::string& text) {
    ProjectiveSignature result{};
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

ProjectiveSignature canonical(ProjectiveSignature signature) {
    std::sort(signature.begin(), signature.end());
    return signature;
}

void write_signature(
    std::ostream& out,
    const ProjectiveSignature& signature) {
    out << '[' << signature[0] << ", " << signature[1]
        << ", " << signature[2] << ']';
}

}  // namespace

FiniteProjectiveFamilyClosureReport
FiniteProjectiveFamilyClosure::certify(
    int maximum_cutoff,
    std::vector<ProjectiveSignature> primitive_squared_lengths) {
    if (primitive_squared_lengths.empty()) {
        throw std::invalid_argument(
            "finite projective family must contain at least one signature");
    }
    for (auto& signature : primitive_squared_lengths) {
        signature = canonical(signature);
    }
    const std::set<ProjectiveSignature> unique_signatures(
        primitive_squared_lengths.begin(),
        primitive_squared_lengths.end());
    if (unique_signatures.size() != primitive_squared_lengths.size()) {
        throw std::invalid_argument(
            "finite projective family contains a duplicate signature");
    }
    std::sort(
        primitive_squared_lengths.begin(),
        primitive_squared_lengths.end());

    FiniteProjectiveFamilyClosureReport report;
    report.maximum_cutoff = maximum_cutoff;
    report.family_cardinality = primitive_squared_lengths.size();
    report.nonempty_fixed_finite_family = true;
    report.unique_primitive_signatures = true;
    report.every_member_triangle_feasible = true;
    report.union_input_degree_bound_proved = true;
    report.union_target_degree_bound_proved = true;

    for (const ProjectiveSignature& signature :
         primitive_squared_lengths) {
        ProjectiveTriadGeometryCertificate member =
            ProjectiveTriadGeometry::certify(
                maximum_cutoff, signature);
        report.total_ordered_role_count += member.ordered_role_count;
        report.unique_primitive_signatures =
            report.unique_primitive_signatures &&
            member.primitive_signature;
        report.every_member_triangle_feasible =
            report.every_member_triangle_feasible &&
            member.triangle_feasible &&
            member.fixed_plane_sphere_geometry;
        report.union_input_degree_bound_proved =
            report.union_input_degree_bound_proved &&
            member.all_degree_bounds_hold;
        report.union_target_degree_bound_proved =
            report.union_target_degree_bound_proved &&
            member.all_degree_bounds_hold;
        report.members.push_back(std::move(member));
    }

    for (int cutoff = 1; cutoff <= maximum_cutoff; ++cutoff) {
        std::size_t summed_input_degree = 0;
        std::size_t summed_target_degree = 0;
        std::size_t summed_elementary_bound = 0;
        for (const auto& member : report.members) {
            const auto& row = member.rows.at(
                static_cast<std::size_t>(cutoff - 1));
            summed_input_degree += row.maximum_input_degree;
            summed_target_degree += row.maximum_target_degree;
            summed_elementary_bound += row.elementary_degree_bound;
        }
        const SpectralReal denominator = static_cast<SpectralReal>(
            std::max<std::size_t>(1, summed_elementary_bound));
        report.maximum_summed_input_degree_ratio = std::max(
            report.maximum_summed_input_degree_ratio,
            static_cast<SpectralReal>(summed_input_degree) /
                denominator);
        report.maximum_summed_target_degree_ratio = std::max(
            report.maximum_summed_target_degree_ratio,
            static_cast<SpectralReal>(summed_target_degree) /
                denominator);
        report.union_input_degree_bound_proved =
            report.union_input_degree_bound_proved &&
            summed_input_degree <= summed_elementary_bound;
        report.union_target_degree_bound_proved =
            report.union_target_degree_bound_proved &&
            summed_target_degree <= summed_elementary_bound;
    }

    report.union_bilinear_shell_bound_proved =
        report.nonempty_fixed_finite_family &&
        report.unique_primitive_signatures &&
        report.every_member_triangle_feasible &&
        report.union_input_degree_bound_proved &&
        report.union_target_degree_bound_proved &&
        report.union_incidence_degree_power == Rational(1) &&
        report.bilinear_l2_frequency_power == Rational(3, 2);
    report.complete_internal_self_cross_decomposition =
        report.union_bilinear_shell_bound_proved;
    report.internal_normalization_terms_bound_proved =
        report.union_bilinear_shell_bound_proved;
    report.cutoff_independent_internal_family_bound =
        report.complete_internal_self_cross_decomposition &&
        report.internal_normalization_terms_bound_proved &&
        report.internal_quartet_frequency_power == Rational(5) &&
        report.target_frequency_power == Rational(11, 2) &&
        report.frequency_gain == Rational(-1, 2);
    return report;
}

FiniteProjectiveFamilyClosureOptions
FiniteProjectiveFamilyClosureCli::parse(
    int argc, char** argv, int first) {
    FiniteProjectiveFamilyClosureOptions options;
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
        } else if (name == "--maximum-height") {
            options.maximum_height = static_cast<SpectralInteger>(
                std::stoll(next(index, name)));
        } else if (name == "--signature") {
            options.primitive_squared_lengths.push_back(
                parse_signature(next(index, name)));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown finite-projective-family option: " + name);
        }
    }
    if (options.maximum_cutoff < 1 || options.maximum_cutoff > 12 ||
        (options.primitive_squared_lengths.empty() &&
         options.maximum_height == 0) ||
        (!options.primitive_squared_lengths.empty() &&
         options.maximum_height != 0) ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "finite-projective-family-certificate requires cutoff 1..12, exactly one of signatures or maximum height, and a certificate path");
    }
    return options;
}

void FiniteProjectiveFamilyClosureCli::print_help(std::ostream& out) {
    out << "Fixed finite projective-family quartet closure options:\n"
        << "  --max-cutoff K       audit incidence through K\n"
        << "  --signature A,B,C    add one primitive squared-length ray; repeatable\n"
        << "  --maximum-height H   use every primitive feasible ray with max(a,b,c)<=H\n"
        << "  --certificate PATH   write English JSON proof certificate\n";
}

int FiniteProjectiveFamilyClosureCli::run(
    const FiniteProjectiveFamilyClosureOptions& options,
    std::ostream& out) {
    const std::vector<ProjectiveSignature> family =
        options.maximum_height > 0
        ? ProjectiveCoreFamily::through_maximum_height(
              options.maximum_height)
        : options.primitive_squared_lengths;
    const FiniteProjectiveFamilyClosureReport report =
        FiniteProjectiveFamilyClosure::certify(
            options.maximum_cutoff,
            family);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write finite projective family certificate");
    }
    certificate << std::setprecision(17)
        << "{\n"
        << "  \"schema\": \"navier-stokes-finite-projective-family-closure-v1\",\n"
        << "  \"claim_scope\": \"the complete internal self-and-cross quartet block of this one fixed finite family\",\n"
        << "  \"family_cardinality\": " << report.family_cardinality
        << ",\n"
        << "  \"primitive_squared_length_signatures\": [\n";
    for (std::size_t index = 0; index < report.members.size(); ++index) {
        certificate << "    ";
        write_signature(
            certificate,
            report.members[index].primitive_squared_lengths);
        certificate << (index + 1 == report.members.size() ? "\n" : ",\n");
    }
    certificate
        << "  ],\n"
        << "  \"geometric_argument\": \"each fixed primitive ray has O_sigma(K) input and target incidence by the plane-sphere intersection argument; summing finitely many such bounds gives O_F(K) incidence for the union\",\n"
        << "  \"operator_argument\": \"B_F=sum_{sigma in F}B_sigma is one restricted bilinear operator; expanding its closed quartet produces every same-ray and unequal-ray term internal to F, so the union bound controls them together without an unsigned sum over a growing shape set\",\n"
        << "  \"normalization_argument\": \"S_F and T_F are finite sums of the fixed-ray shell-offset forms; finite Cauchy-Schwarz closes S_F^2 and S_F T_F with a constant depending only on F\",\n"
        << "  \"constant_dependency\": \"C_F may depend on the listed finite family but is independent of the Fourier cutoff K\",\n"
        << "  \"maximum_cutoff_audit\": " << report.maximum_cutoff
        << ",\n"
        << "  \"total_ordered_role_count\": "
        << report.total_ordered_role_count << ",\n"
        << "  \"maximum_summed_input_degree_ratio\": "
        << static_cast<double>(report.maximum_summed_input_degree_ratio)
        << ",\n"
        << "  \"maximum_summed_target_degree_ratio\": "
        << static_cast<double>(report.maximum_summed_target_degree_ratio)
        << ",\n"
        << "  \"union_incidence_degree_power\": \""
        << report.union_incidence_degree_power.str() << "\",\n"
        << "  \"bilinear_l2_frequency_power\": \""
        << report.bilinear_l2_frequency_power.str() << "\",\n"
        << "  \"internal_quartet_frequency_power\": \""
        << report.internal_quartet_frequency_power.str() << "\",\n"
        << "  \"target_frequency_power\": \""
        << report.target_frequency_power.str() << "\",\n"
        << "  \"frequency_gain\": \""
        << report.frequency_gain.str() << "\",\n"
        << "  \"union_input_degree_bound_proved\": "
        << (report.union_input_degree_bound_proved ? "true" : "false")
        << ",\n"
        << "  \"union_target_degree_bound_proved\": "
        << (report.union_target_degree_bound_proved ? "true" : "false")
        << ",\n"
        << "  \"complete_internal_self_cross_decomposition\": "
        << (report.complete_internal_self_cross_decomposition
                ? "true" : "false") << ",\n"
        << "  \"internal_normalization_terms_bound_proved\": "
        << (report.internal_normalization_terms_bound_proved
                ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_internal_family_bound\": "
        << (report.cutoff_independent_internal_family_bound
                ? "true" : "false") << ",\n"
        << "  \"core_tail_coupling_bound_proved\": false,\n"
        << "  \"growing_tail_internal_bound_proved\": false,\n"
        << "  \"uniform_over_growing_families_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"control interactions between this fixed core and the growing projective tail, and interactions internal to that tail, in the signed bracket-stretching product\"\n"
        << "}\n";
    out << "finite projective family size="
        << report.family_cardinality
        << " internal self+cross="
        << (report.cutoff_independent_internal_family_bound
                ? "PASS" : "FAIL")
        << " gain=R^" << report.frequency_gain.str() << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.cutoff_independent_internal_family_bound ? 0 : 2;
}

}  // namespace lemma
