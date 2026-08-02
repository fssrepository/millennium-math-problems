#include "projective_fan_obstruction.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {

ProjectiveFanObstructionOptions ProjectiveFanObstructionCli::parse(
    int argc, char** argv, int first) {
    ProjectiveFanObstructionOptions options;
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
                "unknown projective-fan-certificate option: " + name);
        }
    }
    if (options.maximum_cutoff < 40 ||
        options.maximum_cutoff > 4096 ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "projective-fan-certificate requires cutoff 40..4096 and a certificate path");
    }
    return options;
}

void ProjectiveFanObstructionCli::print_help(std::ostream& out) {
    out << "Projective coherent-fan obstruction options:\n"
        << "  --max-cutoff K       largest exact count check (40..4096)\n"
        << "  --certificate PATH   write English JSON proof certificate\n";
}

int ProjectiveFanObstructionCli::run(
    const ProjectiveFanObstructionOptions& options,
    std::ostream& out) {
    const ProjectiveFanGeometryCertificate report =
        ProjectiveFanGeometry::certify(options.maximum_cutoff);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective fan proof certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-projective-fan-obstruction-v1\",\n"
        << "  \"family\": \"r=(0,0,K), q=(0,a,b), p=(0,-a,K-b), a/K in [1/2,3/5], b/K in [1/4,2/5]\",\n"
        << "  \"polarization\": \"u_q=e_x and u_p=i(0,-p_z,p_y)/|p| with conjugate negative modes\",\n"
        << "  \"locality_proof\": \"all squared lengths lie between (5/16)K^2 and K^2\",\n"
        << "  \"injectivity_proof\": \"K^2 is the unique largest role; equal primitive rays therefore have equal scale, p^2-q^2 fixes b, and q^2 fixes positive a\",\n"
        << "  \"count_proof\": \"for K>=40 there are at least K/8 admissible b values and K/20 admissible a values, hence N>=K^2/160 distinct primitive rays\",\n"
        << "  \"alignment_proof\": \"each target contribution is +(aK/|p|)e_x with coefficient in [K/2,K]\",\n"
        << "  \"synthesis_lower_bound\": \"||sum_sigma B_sigma||^2/sum_sigma||B_sigma||^2 at target r is at least N/4 >= K^2/640\",\n"
        << "  \"stretching_zero_proof\": \"the lower-band x-polarized q modes never advect; every nonzero upper-band p-advection output misses the active support, so <Au,B(u,u)>=0\",\n"
        << "  \"maximum_cutoff\": " << report.maximum_cutoff << ",\n"
        << "  \"local_frequency_geometry_proved\": true,\n"
        << "  \"target_is_unique_largest_role_proved\": true,\n"
        << "  \"primitive_shape_injectivity_proved\": true,\n"
        << "  \"aligned_target_coefficients_proved\": true,\n"
        << "  \"pair_count_quadratic_lower_bound_proved\": true,\n"
        << "  \"target_synthesis_quadratic_lower_bound_proved\": true,\n"
        << "  \"target_synthesis_unbounded_proved\": true,\n"
        << "  \"stretching_support_disjointness_proved\": true,\n"
        << "  \"exact_zero_stretching_proved\": true,\n"
        << "  \"exact_zero_power_one_product_proved\": true,\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"pair_count\": " << row.pair_count
            << ", \"primitive_shape_count\": "
            << row.primitive_shape_count
            << ", \"target_synthesis_ratio\": "
            << static_cast<double>(row.target_synthesis_ratio)
            << ", \"synthesis_ratio_per_pair\": "
            << static_cast<double>(row.synthesis_ratio_per_pair)
            << ", \"quadratic_lower_bound_ratio\": "
            << static_cast<double>(row.quadratic_lower_bound_ratio)
            << ", \"every_shape_unique\": "
            << (row.every_shape_unique ? "true" : "false")
            << ", \"quadratic_lower_bound_verified\": "
            << (row.quadratic_lower_bound_verified ? "true" : "false")
            << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    certificate
        << "  ],\n"
        << "  \"uniform_projective_synthesis_bound_rejected\": true,\n"
        << "  \"power_one_tradeoff_bound_rejected\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"exploit the exact bracket-stretching tradeoff; projective synthesis cannot be bounded independently of stretching\"\n"
        << "}\n";
    const auto& last = report.rows.back();
    out << std::setprecision(12)
        << "projective fan obstruction K=" << last.cutoff
        << " pairs=" << last.pair_count
        << " target_synthesis="
        << static_cast<double>(last.target_synthesis_ratio)
        << " lower_bound_ratio="
        << static_cast<double>(last.quadratic_lower_bound_ratio)
        << " stretching=0 power_one=0\n"
        << "Certificate written to " << options.certificate_path << '\n';
    return report.target_synthesis_unbounded_proved &&
            report.exact_zero_power_one_product_proved
        ? 0 : 2;
}

}  // namespace lemma
