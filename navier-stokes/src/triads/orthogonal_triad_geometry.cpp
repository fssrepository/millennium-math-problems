#include "orthogonal_triad_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

SpectralInteger integer_dot(WaveVector left, WaveVector right) {
    return static_cast<SpectralInteger>(left.x) * right.x +
        static_cast<SpectralInteger>(left.y) * right.y +
        static_cast<SpectralInteger>(left.z) * right.z;
}

bool inside_cutoff(WaveVector wave, int cutoff) {
    return std::abs(wave.x) <= cutoff &&
        std::abs(wave.y) <= cutoff &&
        std::abs(wave.z) <= cutoff;
}

}  // namespace

OrthogonalTriadGeometryCertificate OrthogonalTriadGeometry::certify(
    int maximum_cutoff) {
    if (maximum_cutoff < 1 || maximum_cutoff > 12) {
        throw std::invalid_argument(
            "orthogonal triad cutoff must be between 1 and 12");
    }
    OrthogonalTriadGeometryCertificate certificate;
    certificate.maximum_cutoff = maximum_cutoff;
    for (int cutoff = 1; cutoff <= maximum_cutoff; ++cutoff) {
        std::map<SpectralInteger, std::vector<WaveVector>> length_groups;
        for (int z = -cutoff; z <= cutoff; ++z) {
            for (int y = -cutoff; y <= cutoff; ++y) {
                for (int x = -cutoff; x <= cutoff; ++x) {
                    if (x == 0 && y == 0 && z == 0) {
                        continue;
                    }
                    const WaveVector wave{x, y, z};
                    length_groups[norm_squared(wave)].push_back(wave);
                }
            }
        }
        OrthogonalTriadCountRow row;
        row.cutoff = cutoff;
        row.elementary_degree_bound =
            2U * static_cast<std::size_t>(2 * cutoff + 1);
        std::map<WaveVector, std::size_t> target_degree;
        for (const auto& [length, waves] : length_groups) {
            static_cast<void>(length);
            for (const WaveVector p : waves) {
                std::size_t input_degree = 0;
                for (const WaveVector q : waves) {
                    if (integer_dot(p, q) != 0) {
                        continue;
                    }
                    const WaveVector target = p + q;
                    if (!inside_cutoff(target, cutoff)) {
                        continue;
                    }
                    ++input_degree;
                    ++target_degree[target];
                    ++row.ordered_pairs;
                }
                row.maximum_input_degree = std::max(
                    row.maximum_input_degree, input_degree);
            }
        }
        for (const auto& [target, degree] : target_degree) {
            static_cast<void>(target);
            row.maximum_target_degree = std::max(
                row.maximum_target_degree, degree);
        }
        certificate.maximum_input_degree_ratio = std::max(
            certificate.maximum_input_degree_ratio,
            static_cast<SpectralReal>(row.maximum_input_degree) /
                static_cast<SpectralReal>(row.elementary_degree_bound));
        certificate.maximum_target_degree_ratio = std::max(
            certificate.maximum_target_degree_ratio,
            static_cast<SpectralReal>(row.maximum_target_degree) /
                static_cast<SpectralReal>(row.elementary_degree_bound));
        certificate.all_degree_bounds_hold =
            certificate.all_degree_bounds_hold &&
            row.maximum_input_degree <= row.elementary_degree_bound &&
            row.maximum_target_degree <= row.elementary_degree_bound;
        certificate.rows.push_back(row);
    }
    return certificate;
}

OrthogonalTriadClosure OrthogonalTriadGeometry::analyze_closure() {
    OrthogonalTriadClosure result;
    result.transfer_frequency_power =
        result.derivative_weight_power +
        result.cauchy_degree_power * result.interaction_degree_power;
    result.generic_local_transfer_frequency_power =
        result.derivative_weight_power +
        result.cauchy_degree_power * result.generic_local_degree_power;
    result.critical_transfer_frequency_power =
        result.derivative_weight_power +
        result.cauchy_degree_power * result.critical_degree_power;
    result.transfer_to_viscosity_frequency_power =
        result.transfer_frequency_power - result.viscous_frequency_power;
    result.transfer_to_viscosity_energy_power =
        result.transfer_shell_energy_power -
        result.viscous_shell_energy_power;
    result.high_frequency_absorbable_from_energy =
        result.transfer_to_viscosity_frequency_power < Rational(0) &&
        result.transfer_to_viscosity_energy_power == Rational(1, 2);
    result.orthogonal_degree_is_subcritical =
        result.interaction_degree_power < result.critical_degree_power;
    result.generic_local_degree_is_supercritical =
        result.critical_degree_power < result.generic_local_degree_power;
    return result;
}

OrthogonalTriadCliOptions OrthogonalTriadCli::parse(
    int argc, char** argv, int first) {
    OrthogonalTriadCliOptions options;
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
                "unknown orthogonal-triad option: " + name);
        }
    }
    return options;
}

void OrthogonalTriadCli::print_help(std::ostream& out) {
    out << "Orthogonal-triad certificate options:\n"
        << "  --max-cutoff K       enumerate cube cutoffs 1 through K\n"
        << "  --certificate PATH   write lattice/closure JSON\n";
}

int OrthogonalTriadCli::run(
    const OrthogonalTriadCliOptions& options, std::ostream& out) {
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "orthogonal-triad-certificate requires --certificate");
    }
    const OrthogonalTriadGeometryCertificate geometry =
        OrthogonalTriadGeometry::certify(options.maximum_cutoff);
    const OrthogonalTriadClosure closure =
        OrthogonalTriadGeometry::analyze_closure();
    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write orthogonal-triad certificate: " +
            options.certificate_path);
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-orthogonal-triad-v1\",\n"
        << "  \"maximum_cutoff\": " << geometry.maximum_cutoff << ",\n"
        << "  \"elementary_degree_bound\": \"2(2K+1)\",\n"
        << "  \"maximum_input_degree_ratio\": "
        << static_cast<double>(geometry.maximum_input_degree_ratio)
        << ",\n  \"maximum_target_degree_ratio\": "
        << static_cast<double>(geometry.maximum_target_degree_ratio)
        << ",\n  \"all_degree_bounds_hold\": "
        << (geometry.all_degree_bounds_hold ? "true" : "false") << ",\n"
        << "  \"closure\": {\"transfer_frequency_power\": \""
        << closure.transfer_frequency_power.str()
        << "\", \"orthogonal_degree_power\": \""
        << closure.interaction_degree_power.str()
        << "\", \"generic_local_degree_power\": \""
        << closure.generic_local_degree_power.str()
        << "\", \"critical_degree_power\": \""
        << closure.critical_degree_power.str()
        << "\", \"generic_local_transfer_frequency_power\": \""
        << closure.generic_local_transfer_frequency_power.str()
        << "\", \"critical_transfer_frequency_power\": \""
        << closure.critical_transfer_frequency_power.str()
        << "\", \"transfer_shell_energy_power\": \""
        << closure.transfer_shell_energy_power.str()
        << "\", \"viscous_frequency_power\": \""
        << closure.viscous_frequency_power.str()
        << "\", \"transfer_to_viscosity_frequency_power\": \""
        << closure.transfer_to_viscosity_frequency_power.str()
        << "\", \"transfer_to_viscosity_energy_power\": \""
        << closure.transfer_to_viscosity_energy_power.str()
        << "\", \"high_frequency_absorbable_from_energy\": "
        << (closure.high_frequency_absorbable_from_energy
                ? "true"
                : "false")
        << ", \"orthogonal_degree_is_subcritical\": "
        << (closure.orthogonal_degree_is_subcritical ? "true" : "false")
        << ", \"generic_local_degree_is_supercritical\": "
        << (closure.generic_local_degree_is_supercritical ? "true" : "false")
        << "},\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < geometry.rows.size(); ++index) {
        const OrthogonalTriadCountRow& row = geometry.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"ordered_pairs\": " << row.ordered_pairs
            << ", \"maximum_input_degree\": "
            << row.maximum_input_degree
            << ", \"maximum_target_degree\": "
            << row.maximum_target_degree
            << ", \"elementary_degree_bound\": "
            << row.elementary_degree_bound << "}"
            << (index + 1 == geometry.rows.size() ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << "orthogonal triad geometry K=1.."
        << options.maximum_cutoff << '\n'
        << "cutoff,ordered_pairs,max_input_degree,max_target_degree,bound\n";
    for (const OrthogonalTriadCountRow& row : geometry.rows) {
        out << row.cutoff << ',' << row.ordered_pairs << ','
            << row.maximum_input_degree << ',' << row.maximum_target_degree
            << ',' << row.elementary_degree_bound << '\n';
    }
    out << "closure transfer K^" << closure.transfer_frequency_power.str()
        << " E^" << closure.transfer_shell_energy_power.str()
        << ", transfer/viscosity K^"
        << closure.transfer_to_viscosity_frequency_power.str()
        << " E^" << closure.transfer_to_viscosity_energy_power.str()
        << ", degree bound="
        << (geometry.all_degree_bounds_hold ? "PASS" : "FAIL") << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return geometry.all_degree_bounds_hold &&
                   closure.high_frequency_absorbable_from_energy
        ? 0
        : 2;
}

}  // namespace lemma
