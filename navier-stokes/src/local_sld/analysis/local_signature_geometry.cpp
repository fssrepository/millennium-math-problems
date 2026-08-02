#include "local_signature_geometry.hpp"

#include "triad_partition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lemma {
namespace {

using OrderedSignature = std::array<SpectralInteger, 3>;
using DegreeKey = std::pair<OrderedSignature, WaveVector>;

SignatureFamilyClosureRow closure_row(Rational signature_count_power) {
    SignatureFamilyClosureRow row;
    row.signature_count_power = signature_count_power;
    row.union_degree_power = Rational(1) + signature_count_power;
    row.transfer_frequency_power = Rational(3) +
        Rational(1, 2) * row.union_degree_power;
    row.transfer_to_viscosity_frequency_power =
        row.transfer_frequency_power - Rational(4);
    row.energy_level_high_frequency_absorption =
        row.transfer_to_viscosity_frequency_power < Rational(0);
    return row;
}

bool inside_cutoff(WaveVector wave, int cutoff) {
    return std::abs(wave.x) <= cutoff &&
        std::abs(wave.y) <= cutoff &&
        std::abs(wave.z) <= cutoff;
}

}  // namespace

LocalSignatureGeometryCertificate LocalSignatureGeometry::certify(
    int maximum_cutoff) {
    if (maximum_cutoff < 1 || maximum_cutoff > 8) {
        throw std::invalid_argument(
            "local signature cutoff must be between 1 and 8");
    }
    LocalSignatureGeometryCertificate certificate;
    certificate.maximum_cutoff = maximum_cutoff;
    for (int cutoff = 1; cutoff <= maximum_cutoff; ++cutoff) {
        std::vector<WaveVector> waves;
        for (int z = -cutoff; z <= cutoff; ++z) {
            for (int y = -cutoff; y <= cutoff; ++y) {
                for (int x = -cutoff; x <= cutoff; ++x) {
                    if (x != 0 || y != 0 || z != 0) {
                        waves.push_back({x, y, z});
                    }
                }
            }
        }
        std::map<DegreeKey, std::size_t> input_degrees;
        std::map<DegreeKey, std::size_t> target_degrees;
        std::map<OrderedSignature, bool> signatures;
        for (const WaveVector p : waves) {
            for (const WaveVector q : waves) {
                const WaveVector target = p + q;
                if (!inside_cutoff(target, cutoff) ||
                    (target.x == 0 && target.y == 0 && target.z == 0) ||
                    !TriadPartitioner::is_local(p, q, target)) {
                    continue;
                }
                const OrderedSignature signature{
                    norm_squared(p), norm_squared(q), norm_squared(target)};
                ++input_degrees[{signature, p}];
                ++target_degrees[{signature, target}];
                signatures[signature] = true;
            }
        }
        LocalSignatureCountRow row;
        row.cutoff = cutoff;
        row.ordered_signatures = signatures.size();
        row.ordered_signature_degree_bound =
            2U * static_cast<std::size_t>(2 * cutoff + 1);
        for (const auto& [key, degree] : input_degrees) {
            static_cast<void>(key);
            row.maximum_input_degree = std::max(
                row.maximum_input_degree, degree);
        }
        for (const auto& [key, degree] : target_degrees) {
            static_cast<void>(key);
            row.maximum_target_degree = std::max(
                row.maximum_target_degree, degree);
        }
        certificate.maximum_input_degree_ratio = std::max(
            certificate.maximum_input_degree_ratio,
            static_cast<SpectralReal>(row.maximum_input_degree) /
                static_cast<SpectralReal>(
                    row.ordered_signature_degree_bound));
        certificate.maximum_target_degree_ratio = std::max(
            certificate.maximum_target_degree_ratio,
            static_cast<SpectralReal>(row.maximum_target_degree) /
                static_cast<SpectralReal>(
                    row.ordered_signature_degree_bound));
        certificate.all_fixed_signature_degree_bounds_hold =
            certificate.all_fixed_signature_degree_bounds_hold &&
            row.maximum_input_degree <=
                row.ordered_signature_degree_bound &&
            row.maximum_target_degree <=
                row.ordered_signature_degree_bound;
        certificate.rows.push_back(row);
    }
    return certificate;
}

SignatureFamilyClosure LocalSignatureGeometry::analyze_closure() {
    SignatureFamilyClosure result;
    result.finite_signature_family = closure_row(Rational(0));
    result.critical_signature_family = closure_row(Rational(1));
    result.dense_signature_family = closure_row(Rational(2));
    result.closing_requires_sublinear_signature_count =
        result.finite_signature_family.energy_level_high_frequency_absorption &&
        !result.critical_signature_family.energy_level_high_frequency_absorption &&
        !result.dense_signature_family.energy_level_high_frequency_absorption;
    return result;
}

LocalSignatureCliOptions LocalSignatureCli::parse(
    int argc, char** argv, int first) {
    LocalSignatureCliOptions options;
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
                "unknown local-signature option: " + name);
        }
    }
    return options;
}

void LocalSignatureCli::print_help(std::ostream& out) {
    out << "Local-signature certificate options:\n"
        << "  --max-cutoff K       enumerate local ordered signatures\n"
        << "  --certificate PATH   write degree/closure JSON\n";
}

int LocalSignatureCli::run(
    const LocalSignatureCliOptions& options, std::ostream& out) {
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "local-signature-certificate requires --certificate");
    }
    const LocalSignatureGeometryCertificate geometry =
        LocalSignatureGeometry::certify(options.maximum_cutoff);
    const SignatureFamilyClosure closure =
        LocalSignatureGeometry::analyze_closure();
    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write local-signature certificate: " +
            options.certificate_path);
    }
    auto write_closure_row = [&](const char* name,
                                 const SignatureFamilyClosureRow& row) {
        certificate << "    \"" << name << "\": {"
            << "\"signature_count_power\": \""
            << row.signature_count_power.str()
            << "\", \"union_degree_power\": \""
            << row.union_degree_power.str()
            << "\", \"transfer_frequency_power\": \""
            << row.transfer_frequency_power.str()
            << "\", \"transfer_to_viscosity_frequency_power\": \""
            << row.transfer_to_viscosity_frequency_power.str()
            << "\", \"energy_level_high_frequency_absorption\": "
            << (row.energy_level_high_frequency_absorption
                    ? "true"
                    : "false") << '}';
    };
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-local-signature-v1\",\n"
        << "  \"maximum_cutoff\": " << geometry.maximum_cutoff << ",\n"
        << "  \"ordered_signature_degree_bound\": \"2(2K+1)\",\n"
        << "  \"square_summed_transfer_bound\": "
           "\"(sum_sigma |V_sigma|^2)^(1/2) <= "
           "C K^(7/2) E_K^(3/2)\",\n"
        << "  \"effective_signature_identity\": "
           "\"sum_sigma |V_sigma| = sqrt(N_eff) "
           "(sum_sigma |V_sigma|^2)^(1/2)\",\n"
        << "  \"signed_amplification_identity\": "
           "\"|sum_sigma V_sigma| = A_sig "
           "(sum_sigma |V_sigma|^2)^(1/2)\",\n"
        << "  \"critical_signed_amplification_power\": \""
        << closure.critical_signed_amplification_power.str()
        << "\",\n"
        << "  \"maximum_input_degree_ratio\": "
        << static_cast<double>(geometry.maximum_input_degree_ratio)
        << ",\n  \"maximum_target_degree_ratio\": "
        << static_cast<double>(geometry.maximum_target_degree_ratio)
        << ",\n  \"all_fixed_signature_degree_bounds_hold\": "
        << (geometry.all_fixed_signature_degree_bounds_hold
                ? "true"
                : "false") << ",\n"
        << "  \"closure\": {\n";
    write_closure_row("finite", closure.finite_signature_family);
    certificate << ",\n";
    write_closure_row("critical", closure.critical_signature_family);
    certificate << ",\n";
    write_closure_row("dense", closure.dense_signature_family);
    certificate << "\n  },\n  \"rows\": [\n";
    for (std::size_t index = 0; index < geometry.rows.size(); ++index) {
        const LocalSignatureCountRow& row = geometry.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"ordered_signatures\": " << row.ordered_signatures
            << ", \"maximum_input_degree\": "
            << row.maximum_input_degree
            << ", \"maximum_target_degree\": "
            << row.maximum_target_degree
            << ", \"ordered_signature_degree_bound\": "
            << row.ordered_signature_degree_bound << "}"
            << (index + 1 == geometry.rows.size() ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << "local signature geometry K=1.." << options.maximum_cutoff
        << "\ncutoff,ordered_signatures,max_input_degree,"
           "max_target_degree,bound\n";
    for (const LocalSignatureCountRow& row : geometry.rows) {
        out << row.cutoff << ',' << row.ordered_signatures << ','
            << row.maximum_input_degree << ',' << row.maximum_target_degree
            << ',' << row.ordered_signature_degree_bound << '\n';
    }
    out << "finite signatures transfer K^"
        << closure.finite_signature_family.transfer_frequency_power.str()
        << ", critical count K^"
        << closure.critical_signature_count_power.str()
        << ", critical signed amplification K^"
        << closure.critical_signed_amplification_power.str()
        << ", square-summed bound="
        << (closure.square_summed_fixed_signature_bound ? "PASS" : "FAIL")
        << ", fixed-signature degree bound="
        << (geometry.all_fixed_signature_degree_bounds_hold
                ? "PASS"
                : "FAIL") << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return geometry.all_fixed_signature_degree_bounds_hold &&
                   closure.closing_requires_sublinear_signature_count
        ? 0
        : 2;
}

}  // namespace lemma
