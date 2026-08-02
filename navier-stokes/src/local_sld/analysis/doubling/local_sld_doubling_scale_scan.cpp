#include "local_sld_doubling_scale_scan.hpp"

#include "local_quartic_closure_objective.hpp"
#include "local_sld_two_scale_state.hpp"
#include "local_sld_signature_block.hpp"
#include "spectral_galerkin.hpp"

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

SpectralReal target_scale(SpectralReal enstrophy,
                          SpectralReal palinstrophy) {
    return std::pow(enstrophy, 1.25L) *
        std::pow(palinstrophy, 0.75L);
}

void write_row(std::ostream& output,
               const LocalSldDoublingScaleRow& row) {
    output << "{\"scale\": " << row.scale
        << ", \"response_angle\": "
        << static_cast<double>(row.response_angle)
        << ", \"high_to_low_energy_ratio\": "
        << static_cast<double>(row.high_to_low_energy_ratio)
        << ", \"enstrophy\": " << static_cast<double>(row.enstrophy)
        << ", \"palinstrophy\": "
        << static_cast<double>(row.palinstrophy)
        << ", \"bracket\": " << static_cast<double>(row.bracket)
        << ", \"signed_target_ratio\": "
        << static_cast<double>(row.signed_target_ratio)
        << ", \"full_local_bracket\": "
        << static_cast<double>(row.full_local_bracket)
        << ", \"full_local_target_ratio\": "
        << static_cast<double>(row.full_local_target_ratio)
        << ", \"remainder_closed_bracket\": "
        << static_cast<double>(row.remainder_closed_bracket)
        << ", \"remainder_closed_target_ratio\": "
        << static_cast<double>(row.remainder_closed_target_ratio)
        << ", \"mixed_bracket\": "
        << static_cast<double>(row.mixed_bracket)
        << ", \"mixed_target_ratio\": "
        << static_cast<double>(row.mixed_target_ratio)
        << ", \"projected\": " << static_cast<double>(row.projected)
        << ", \"projected_target_ratio\": "
        << static_cast<double>(row.projected_target_ratio)
        << ", \"non_projected\": "
        << static_cast<double>(row.non_projected)
        << ", \"non_projected_target_ratio\": "
        << static_cast<double>(row.non_projected_target_ratio)
        << ", \"cancellation_fraction\": "
        << static_cast<double>(row.cancellation_fraction)
        << ", \"square_identity_error\": "
        << static_cast<double>(row.square_identity_error) << '}';
}

void write_certificate(
    const LocalSldDoublingScaleScanReport& report,
    const LocalSldDoublingScaleScanOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write doubling scale-scan certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-doubling-scale-scan-v1\",\n"
        << "  \"definition\": \"two-scale signed screen of the complete closed (m,m,2m) K+G bracket\",\n"
        << "  \"minimum_scale\": " << report.minimum_scale << ",\n"
        << "  \"maximum_scale\": " << report.maximum_scale << ",\n"
        << "  \"minimum_angle\": "
        << static_cast<double>(options.minimum_angle) << ",\n"
        << "  \"maximum_angle\": "
        << static_cast<double>(options.maximum_angle) << ",\n"
        << "  \"angle_count\": " << report.angle_count << ",\n"
        << "  \"energy_decay_power\": "
        << static_cast<double>(report.energy_decay_power) << ",\n"
        << "  \"maximum_signed_row\": ";
    write_row(output, report.maximum_signed_row);
    output << ",\n  \"maximum_absolute_row\": ";
    write_row(output, report.maximum_absolute_row);
    output << ",\n  \"maximum_non_projected_row\": ";
    write_row(output, report.maximum_non_projected_row);
    output << ",\n  \"maximum_remainder_row\": ";
    write_row(output, report.maximum_remainder_row);
    output << ",\n  \"maximum_mixed_row\": ";
    write_row(output, report.maximum_mixed_row);
    output << ",\n  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        output << "    ";
        write_row(output, report.rows[index]);
        output << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_square_identity_verified\": "
        << (report.every_square_identity_verified ? "true" : "false")
        << ",\n"
        << "  \"every_block_decomposition_verified\": "
        << (report.every_block_decomposition_verified
            ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_bound_proved\": false,\n"
        << "  \"finite_scan_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldDoublingScaleScanReport LocalSldDoublingScaleScan::analyze(
    const LocalSldDoublingScaleScanOptions& options) {
    SpectralGalerkin galerkin;
    galerkin.configure("direct", options.threads);
    const SpectralDynamics dynamics(galerkin);
    const TriadSelection selection =
        TriadSelection::local_equal_low_doubling();
    std::vector<SpectralReal> angles;
    angles.reserve(static_cast<std::size_t>(options.angle_count));
    for (int index = 0; index < options.angle_count; ++index) {
        const SpectralReal fraction = options.angle_count == 1
            ? 0.0L
            : static_cast<SpectralReal>(index) /
                  static_cast<SpectralReal>(options.angle_count - 1);
        angles.push_back(options.minimum_angle + fraction *
            (options.maximum_angle - options.minimum_angle));
    }

    LocalSldDoublingScaleScanReport report;
    report.minimum_scale = options.minimum_scale;
    report.maximum_scale = options.maximum_scale;
    report.angle_count = options.angle_count;
    report.energy_decay_power = options.energy_decay_power;
    report.every_square_identity_verified = true;
    report.every_block_decomposition_verified = true;
    bool first = true;
    for (int scale = options.minimum_scale;
         scale <= options.maximum_scale; ++scale) {
        const SpectralReal energy_ratio = std::pow(
            static_cast<SpectralReal>(scale),
            -options.energy_decay_power);
        const std::vector<SpectralState> states =
            LocalSldTwoScaleState::cyclic_response_mixtures(
                dynamics, scale, energy_ratio, angles);
        for (std::size_t index = 0; index < states.size(); ++index) {
            const LocalSldSignatureBlockReport block =
                LocalSldSignatureBlock::analyze(
                    dynamics, states[index], {1, 1, 2}, true);
            const LocalQuarticClosureObjectiveValue& full =
                block.dominant;
            const LocalSldProjectedSquareReport projected =
                LocalSldProjectedSquare::evaluate(
                    dynamics, states[index], selection);
            LocalSldDoublingScaleRow row;
            row.scale = scale;
            row.response_angle = angles[index];
            row.high_to_low_energy_ratio = energy_ratio;
            row.enstrophy = full.enstrophy;
            row.palinstrophy = full.palinstrophy;
            row.bracket = full.signed_two_entry_bracket;
            row.full_local_bracket = block.full_bracket;
            row.remainder_closed_bracket =
                block.remainder_closed_bracket;
            row.mixed_bracket = block.cross_bracket;
            row.projected = projected.expanded_total;
            row.non_projected = row.bracket - row.projected;
            const SpectralReal denominator = target_scale(
                row.enstrophy, row.palinstrophy);
            if (denominator > 0.0L) {
                row.signed_target_ratio = row.bracket / denominator;
                row.full_local_target_ratio =
                    row.full_local_bracket / denominator;
                row.remainder_closed_target_ratio =
                    row.remainder_closed_bracket / denominator;
                row.mixed_target_ratio =
                    row.mixed_bracket / denominator;
                row.projected_target_ratio =
                    row.projected / denominator;
                row.non_projected_target_ratio =
                    row.non_projected / denominator;
            }
            const SpectralReal absolute_parts =
                std::abs(row.projected) + std::abs(row.non_projected);
            if (absolute_parts > 0.0L) {
                row.cancellation_fraction = 1.0L -
                    std::abs(row.bracket) / absolute_parts;
            }
            row.square_identity_error =
                projected.completion_relative_error;
            report.every_square_identity_verified =
                report.every_square_identity_verified &&
                projected.identity_verified;
            report.every_block_decomposition_verified =
                report.every_block_decomposition_verified &&
                block.exact_decomposition;
            if (first || row.signed_target_ratio >
                             report.maximum_signed_row.signed_target_ratio) {
                report.maximum_signed_row = row;
            }
            if (first || std::abs(row.signed_target_ratio) >
                             std::abs(report.maximum_absolute_row
                                          .signed_target_ratio)) {
                report.maximum_absolute_row = row;
            }
            if (first || std::abs(row.non_projected_target_ratio) >
                             std::abs(report.maximum_non_projected_row
                                          .non_projected_target_ratio)) {
                report.maximum_non_projected_row = row;
            }
            if (first || std::abs(row.remainder_closed_target_ratio) >
                             std::abs(report.maximum_remainder_row
                                          .remainder_closed_target_ratio)) {
                report.maximum_remainder_row = row;
            }
            if (first || std::abs(row.mixed_target_ratio) >
                             std::abs(report.maximum_mixed_row
                                          .mixed_target_ratio)) {
                report.maximum_mixed_row = row;
            }
            first = false;
            report.rows.push_back(row);
        }
    }
    return report;
}

LocalSldDoublingScaleScanOptions LocalSldDoublingScaleScan::parse(
    int argc, char** argv, int first) {
    LocalSldDoublingScaleScanOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--min-scale") {
            options.minimum_scale = std::stoi(next(index, name));
        } else if (name == "--max-scale") {
            options.maximum_scale = std::stoi(next(index, name));
        } else if (name == "--angle-min") {
            options.minimum_angle = std::stold(next(index, name));
        } else if (name == "--angle-max") {
            options.maximum_angle = std::stold(next(index, name));
        } else if (name == "--angle-count") {
            options.angle_count = std::stoi(next(index, name));
        } else if (name == "--energy-decay-power") {
            options.energy_decay_power = std::stold(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown doubling scale-scan option: " + name);
        }
    }
    if (options.minimum_scale < 2 ||
        options.maximum_scale < options.minimum_scale ||
        options.maximum_scale > 12 || options.angle_count < 1 ||
        options.angle_count > 1001 ||
        !std::isfinite(options.minimum_angle) ||
        !std::isfinite(options.maximum_angle) ||
        options.maximum_angle < options.minimum_angle ||
        !(options.energy_decay_power > 0.0L) ||
        !std::isfinite(options.energy_decay_power) ||
        options.threads < 1 || options.threads > 256 ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "doubling scale scan requires scales 2..12, finite angles, positive decay power, and a certificate");
    }
    return options;
}

void LocalSldDoublingScaleScan::print_help(std::ostream& out) {
    out << "Local SLD doubling two-scale scan options:\n"
        << "  --min-scale L        minimum integer dilation (default 2)\n"
        << "  --max-scale L        maximum integer dilation (default 12)\n"
        << "  --angle-min X        minimum axis/response angle\n"
        << "  --angle-max X        maximum axis/response angle\n"
        << "  --angle-count N      inclusive uniform angle count\n"
        << "  --energy-decay-power X  high energy L^(-X), default 2.75\n"
        << "  --threads N          direct-kernel workers\n"
        << "  --certificate PATH   write English JSON scan\n";
}

int LocalSldDoublingScaleScan::run(
    const LocalSldDoublingScaleScanOptions& options,
    std::ostream& out) {
    const LocalSldDoublingScaleScanReport report = analyze(options);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "doubling scale scan rows=" << report.rows.size()
        << " max_signed="
        << static_cast<double>(
               report.maximum_signed_row.signed_target_ratio)
        << " at L=" << report.maximum_signed_row.scale
        << " angle="
        << static_cast<double>(report.maximum_signed_row.response_angle)
        << " max_abs="
        << static_cast<double>(
               report.maximum_absolute_row.signed_target_ratio)
        << " square_identities="
        << (report.every_square_identity_verified ? "PASS" : "FAIL")
        << " block_decompositions="
        << (report.every_block_decomposition_verified ? "PASS" : "FAIL")
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_square_identity_verified &&
            report.every_block_decomposition_verified
        ? 0 : 2;
}

}  // namespace lemma
