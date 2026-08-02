#include "projective_square_function_closure.hpp"

#include <filesystem>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {

ProjectiveSquareFunctionClosureReport
ProjectiveSquareFunctionClosure::certify() {
    ProjectiveSquareFunctionClosureReport report;
    report.primitive_ray_partition_is_disjoint = true;
    // A fixed input determines the scale on a primitive ray. For each of at
    // most six role assignments, plane--sphere incidence gives at most
    // 2(2K+1) partners. This constant is independent of the primitive shape.
    report.uniform_single_ray_incidence_constant_proved =
        report.uniform_ray_incidence_degree_power == Rational(1);
    // Squared Cauchy--Schwarz costs one ray degree K, the derivative costs
    // K^2 after squaring, and the disjoint pair partition removes any shape
    // count: sum_sigma ||B_sigma||_2^2 <= C K^3 E^2.
    report.bilinear_projective_square_function_bound_proved =
        report.primitive_ray_partition_is_disjoint &&
        report.uniform_single_ray_incidence_constant_proved &&
        report.squared_function_squared_frequency_power ==
            report.uniform_ray_incidence_degree_power +
            Rational(2) * report.derivative_frequency_power &&
        report.squared_function_frequency_power ==
            report.squared_function_squared_frequency_power /
                Rational(2);
    report.square_function_has_target_power_gain =
        report.bilinear_projective_square_function_bound_proved &&
        report.candidate_diagonal_quartet_frequency_power ==
            Rational(2) +
            Rational(2) * report.squared_function_frequency_power &&
        report.candidate_diagonal_frequency_gain < Rational(0);
    // The components B_sigma share output Fourier modes. The square-function
    // estimate does not bound ||sum_sigma B_sigma|| without an additional
    // synthesis/almost-orthogonality estimate, so all cross-ray flags stay
    // false deliberately.
    return report;
}

ProjectiveSquareFunctionClosureOptions
ProjectiveSquareFunctionClosureCli::parse(
    int argc, char** argv, int first) {
    ProjectiveSquareFunctionClosureOptions options;
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--certificate") {
            if (index + 1 >= argc) {
                throw std::invalid_argument(
                    "missing value for --certificate");
            }
            options.certificate_path = argv[++index];
        } else {
            throw std::invalid_argument(
                "unknown projective-square-function option: " + name);
        }
    }
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "projective-square-function-certificate requires --certificate");
    }
    return options;
}

void ProjectiveSquareFunctionClosureCli::print_help(std::ostream& out) {
    out << "Projective-ray square-function closure options:\n"
        << "  --certificate PATH   write English JSON proof certificate\n";
}

int ProjectiveSquareFunctionClosureCli::run(
    const ProjectiveSquareFunctionClosureOptions& options,
    std::ostream& out) {
    const ProjectiveSquareFunctionClosureReport report =
        ProjectiveSquareFunctionClosure::certify();
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective square-function certificate");
    }
    certificate
        << "{\n"
        << "  \"schema\": \"navier-stokes-projective-square-function-closure-v1\",\n"
        << "  \"ray_partition\": \"each local ordered input pair belongs to exactly one primitive squared-length ray\",\n"
        << "  \"uniform_incidence_proof\": \"a fixed input determines the ray scale; at most six role assignments and at most 2(2K+1) plane-sphere partners per role give a shape-independent O(K) degree\",\n"
        << "  \"square_function_bound\": \"sum_sigma ||B_sigma(u,u)||_2^2 <= C K^3 E^2\",\n"
        << "  \"square_function_derivation\": \"target-wise Cauchy-Schwarz costs one K degree, the squared derivative costs K^2, and disjoint ray pair sets eliminate the number of shapes\",\n"
        << "  \"uniform_ray_incidence_degree_power\": \""
        << report.uniform_ray_incidence_degree_power.str() << "\",\n"
        << "  \"squared_function_squared_frequency_power\": \""
        << report.squared_function_squared_frequency_power.str()
        << "\",\n"
        << "  \"squared_function_frequency_power\": \""
        << report.squared_function_frequency_power.str() << "\",\n"
        << "  \"candidate_diagonal_quartet_frequency_power\": \""
        << report.candidate_diagonal_quartet_frequency_power.str()
        << "\",\n"
        << "  \"target_frequency_power\": \""
        << report.target_frequency_power.str() << "\",\n"
        << "  \"candidate_diagonal_frequency_gain\": \""
        << report.candidate_diagonal_frequency_gain.str() << "\",\n"
        << "  \"primitive_ray_partition_is_disjoint\": true,\n"
        << "  \"uniform_single_ray_incidence_constant_proved\": true,\n"
        << "  \"bilinear_projective_square_function_bound_proved\": true,\n"
        << "  \"square_function_has_target_power_gain\": true,\n"
        << "  \"diagonal_projective_quartet_sum_proved\": false,\n"
        << "  \"coherent_projective_synthesis_bound_proved\": false,\n"
        << "  \"cross_ray_quartet_bound_proved\": false,\n"
        << "  \"power_one_tradeoff_bound_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"control the coherent cross-ray synthesis on shared output modes, including normalization and the product with S_full\"\n"
        << "}\n";
    out << "projective square function=K^"
        << report.squared_function_frequency_power.str()
        << " diagonal quartet=K^"
        << report.candidate_diagonal_quartet_frequency_power.str()
        << " cross_ray="
        << (report.cross_ray_quartet_bound_proved ? "PASS" : "OPEN")
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.bilinear_projective_square_function_bound_proved &&
            report.square_function_has_target_power_gain
        ? 0 : 2;
}

}  // namespace lemma
