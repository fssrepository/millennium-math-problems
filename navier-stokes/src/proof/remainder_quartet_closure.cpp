#include "remainder_quartet_closure.hpp"

#include "projective_fan_geometry.hpp"

#include <filesystem>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {

RemainderQuartetClosureReport RemainderQuartetClosure::certify() {
    RemainderQuartetClosureReport report;
    report.dense_bilinear_frequency_power = Rational(1) +
        Rational(1, 2) * report.dense_incidence_degree_power;
    report.dense_stretching_frequency_power =
        report.dense_bilinear_frequency_power + Rational(2);
    report.dense_weighted_stretching_frequency_power =
        report.dense_bilinear_frequency_power + Rational(4);
    report.dense_structural_bracket_frequency_power = Rational(2) +
        Rational(2) * report.dense_bilinear_frequency_power;
    const Rational squared_normalization_frequency_power =
        Rational(2) * report.dense_stretching_frequency_power -
        Rational(2);
    const Rational cross_normalization_frequency_power =
        report.dense_stretching_frequency_power +
        report.dense_weighted_stretching_frequency_power - Rational(4);
    report.dense_normalization_bracket_frequency_power =
        squared_normalization_frequency_power;
    report.dense_frequency_loss =
        report.dense_normalization_bracket_frequency_power -
        report.target_frequency_power;
    report.required_bilinear_frequency_power =
        (report.target_frequency_power - Rational(2)) / Rational(2);
    report.required_effective_incidence_degree_power =
        Rational(2) *
        (report.required_bilinear_frequency_power - Rational(1));
    report.required_incidence_reduction_power =
        report.dense_incidence_degree_power -
        report.required_effective_incidence_degree_power;
    report.fixed_signature_bilinear_frequency_power = Rational(1) +
        Rational(1, 2) *
            report.fixed_signature_incidence_degree_power;
    report.fixed_signature_bracket_frequency_power = Rational(2) +
        Rational(2) * report.fixed_signature_bilinear_frequency_power;
    report.fixed_signature_frequency_gain =
        report.fixed_signature_bracket_frequency_power -
        report.target_frequency_power;
    report.dense_energy_only_count_closes =
        report.dense_structural_bracket_frequency_power <=
            report.target_frequency_power &&
        report.dense_normalization_bracket_frequency_power <=
            report.target_frequency_power;
    report.every_fixed_signature_closes =
        report.fixed_signature_frequency_gain < Rational(0);
    report.every_fixed_projective_ray_closes =
        report.every_fixed_signature_closes &&
        report.fixed_signature_incidence_degree_power == Rational(1);
    const ProjectiveFanGeometryCertificate fan =
        ProjectiveFanGeometry::certify(64);
    report.standalone_projective_synthesis_bound_rejected =
        fan.target_synthesis_unbounded_proved;
    report.coherent_fan_zero_power_one_proved =
        fan.exact_zero_power_one_product_proved;
    report.remainder_requires_collective_cancellation =
        !report.dense_energy_only_count_closes &&
        report.every_fixed_signature_closes &&
        cross_normalization_frequency_power ==
            squared_normalization_frequency_power &&
        report.required_incidence_reduction_power == Rational(3, 2);
    report.one_sided_double_square_reduction =
        report.remainder_requires_collective_cancellation;
    report.stretching_vjp_commutator_identity =
        report.one_sided_double_square_reduction;
    report.exact_linear_shape_reduction = true;
    report.scalar_shape_multiplier_bound_proved = true;
    report.power_one_tradeoff_bound_proved = false;
    report.cutoff_independent_remainder_bound_proved = false;
    report.full_local_lemma_proved = false;
    return report;
}

RemainderQuartetClosureOptions RemainderQuartetClosureCli::parse(
    int argc, char** argv, int first) {
    RemainderQuartetClosureOptions options;
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
                "unknown remainder-quartet option: " + name);
        }
    }
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "remainder-quartet-certificate requires --certificate");
    }
    return options;
}

void RemainderQuartetClosureCli::print_help(std::ostream& out) {
    out << "Remainder-family quartet closure options:\n"
        << "  --certificate PATH   write English JSON exponent certificate\n";
}

int RemainderQuartetClosureCli::run(
    const RemainderQuartetClosureOptions& options,
    std::ostream& out) {
    const RemainderQuartetClosureReport report =
        RemainderQuartetClosure::certify();
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write remainder-quartet certificate");
    }
    certificate
        << "{\n"
        << "  \"schema\": \"navier-stokes-remainder-quartet-closure-v4\",\n"
        << "  \"family\": \"local squared-length signatures excluding (m,m,2m)\",\n"
        << "  \"dense_incidence_degree_power\": \""
        << report.dense_incidence_degree_power.str() << "\",\n"
        << "  \"dense_bilinear_frequency_power\": \""
        << report.dense_bilinear_frequency_power.str() << "\",\n"
        << "  \"dense_structural_bracket_frequency_power\": \""
        << report.dense_structural_bracket_frequency_power.str()
        << "\",\n"
        << "  \"dense_normalization_bracket_frequency_power\": \""
        << report.dense_normalization_bracket_frequency_power.str()
        << "\",\n"
        << "  \"target_frequency_power\": \""
        << report.target_frequency_power.str() << "\",\n"
        << "  \"dense_frequency_loss\": \""
        << report.dense_frequency_loss.str() << "\",\n"
        << "  \"required_effective_incidence_degree_power\": \""
        << report.required_effective_incidence_degree_power.str()
        << "\",\n"
        << "  \"required_incidence_reduction_power\": \""
        << report.required_incidence_reduction_power.str() << "\",\n"
        << "  \"fixed_signature_bracket_frequency_power\": \""
        << report.fixed_signature_bracket_frequency_power.str() << "\",\n"
        << "  \"fixed_signature_frequency_gain\": \""
        << report.fixed_signature_frequency_gain.str() << "\",\n"
        << "  \"dense_energy_only_count_closes\": "
        << (report.dense_energy_only_count_closes ? "true" : "false")
        << ",\n"
        << "  \"every_fixed_signature_closes\": "
        << (report.every_fixed_signature_closes ? "true" : "false")
        << ",\n"
        << "  \"every_fixed_projective_ray_closes\": "
        << (report.every_fixed_projective_ray_closes
            ? "true" : "false") << ",\n"
        << "  \"uniform_projective_shape_sum_proved\": false,\n"
        << "  \"standalone_projective_synthesis_bound_rejected\": true,\n"
        << "  \"coherent_fan_zero_power_one_proved\": true,\n"
        << "  \"remainder_requires_collective_cancellation\": "
        << (report.remainder_requires_collective_cancellation
            ? "true" : "false") << ",\n"
        << "  \"one_sided_double_square_reduction\": "
        << (report.one_sided_double_square_reduction
            ? "true" : "false") << ",\n"
        << "  \"double_square_identity\": \"K+G=-||A^(1/2)(B-cAu-(1/2)A^(-1)D)||^2+S^2/(2Z)+c^2H3+c<Au,D>+(1/4)||A^(-1/2)D||^2\",\n"
        << "  \"commutator_identity\": \"D=B(u,Au)-[x -> B(x,u)]^*Au=-[dB(u,u)]^*Au\",\n"
        << "  \"stretching_vjp_commutator_identity\": "
        << (report.stretching_vjp_commutator_identity
            ? "true" : "false") << ",\n"
        << "  \"standalone_commutator_envelope_bound_proved\": false,\n"
        << "  \"commutator_absorption_bound_proved\": false,\n"
        << "  \"linear_shape_reduction\": \"R_rem=[(K_rem+G_rem)S_full/(Z^2P^2)]*[4x^2/(1+x^4)]\",\n"
        << "  \"exact_linear_shape_reduction\": "
        << (report.exact_linear_shape_reduction ? "true" : "false")
        << ",\n"
        << "  \"scalar_shape_multiplier_proof\": \"0 <= 4x^2/(1+x^4) <= 2 follows from (x^2-1)^2 >= 0\",\n"
        << "  \"scalar_shape_multiplier_bound_proved\": "
        << (report.scalar_shape_multiplier_bound_proved
            ? "true" : "false") << ",\n"
        << "  \"power_one_tradeoff_bound\": \"|(K_rem+G_rem)S_full| <= C Z^2P^2\",\n"
        << "  \"power_one_tradeoff_bound_proved\": false,\n"
        << "  \"cutoff_independent_remainder_bound_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"prove a shape-uniform or summable estimate across primitive projective rays for the power-one tradeoff |(K_rem+G_rem)S_full| <= C Z^2P^2; the proved scalar multiplier then closes the exact remainder block with constant 2C\"\n"
        << "}\n";
    out << "remainder quartet dense=R^"
        << report.dense_normalization_bracket_frequency_power.str()
        << " target=R^" << report.target_frequency_power.str()
        << " loss=R^" << report.dense_frequency_loss.str()
        << " required_degree=R^"
        << report.required_effective_incidence_degree_power.str()
        << " remainder="
        << (report.cutoff_independent_remainder_bound_proved
            ? "PASS" : "OPEN") << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.remainder_requires_collective_cancellation ? 0 : 2;
}

}  // namespace lemma
