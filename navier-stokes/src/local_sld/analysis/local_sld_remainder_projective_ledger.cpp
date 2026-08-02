#include "local_sld_remainder_projective_ledger.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

struct Accumulator {
    LocalSldRemainderProjectiveEntry entry;
};

SpectralReal relative_error(SpectralReal left, SpectralReal right) {
    return std::abs(left - right) /
        std::max({std::abs(left), std::abs(right), 1e-30L});
}

std::pair<Shape, SpectralInteger> primitive_shape(
    const Shape& signature) {
    const SpectralInteger divisor = std::gcd(
        signature[0], std::gcd(signature[1], signature[2]));
    if (divisor <= 0) {
        throw std::invalid_argument(
            "projective signature requires positive squared lengths");
    }
    return {{signature[0] / divisor,
             signature[1] / divisor,
             signature[2] / divisor}, divisor};
}

void write_json(
    const LocalSldRemainderProjectiveReport& report,
    const LocalSldRemainderProjectiveCliOptions& options,
    std::ostream& output) {
    const std::size_t count = std::min(
        report.shapes.size(), static_cast<std::size_t>(options.top));
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-remainder-projective-ledger-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"signature_count\": " << report.signature_count << ",\n"
        << "  \"projective_shape_count\": "
        << report.projective_shape_count << ",\n"
        << "  \"excludes_signature_123\": "
        << (report.excludes_signature_123 ? "true" : "false")
        << ",\n"
        << "  \"excludes_triple_family\": "
        << (report.excludes_triple_family ? "true" : "false")
        << ",\n"
        << "  \"expected_power_one_total\": "
        << static_cast<double>(report.expected_power_one_total) << ",\n"
        << "  \"reconstructed_power_one_total\": "
        << static_cast<double>(report.reconstructed_power_one_total)
        << ",\n"
        << "  \"reconstruction_error\": "
        << static_cast<double>(report.reconstruction_error) << ",\n"
        << "  \"absolute_projective_sum\": "
        << static_cast<double>(report.absolute_projective_sum) << ",\n"
        << "  \"squared_projective_sum\": "
        << static_cast<double>(report.squared_projective_sum) << ",\n"
        << "  \"effective_projective_shapes\": "
        << static_cast<double>(report.effective_projective_shapes)
        << ",\n"
        << "  \"dominant_projective_fraction\": "
        << static_cast<double>(report.dominant_projective_fraction)
        << ",\n"
        << "  \"signed_projective_alignment\": "
        << static_cast<double>(report.signed_projective_alignment)
        << ",\n"
        << "  \"reported_shape_count\": " << count << ",\n"
        << "  \"top_projective_shapes\": [\n";
    for (std::size_t index = 0; index < count; ++index) {
        const LocalSldRemainderProjectiveEntry& row =
            report.shapes[index];
        output << "    {\"primitive_squared_lengths\": ["
            << row.primitive_squared_lengths[0] << ", "
            << row.primitive_squared_lengths[1] << ", "
            << row.primitive_squared_lengths[2]
            << "], \"scale_count\": " << row.scale_count
            << ", \"interactions\": " << row.interactions
            << ", \"minimum_scale\": " << row.minimum_scale
            << ", \"maximum_scale\": " << row.maximum_scale
            << ", \"bracket_total\": "
            << static_cast<double>(row.bracket_total)
            << ", \"power_one_total\": "
            << static_cast<double>(row.power_one_total)
            << ", \"absolute_scale_sum\": "
            << static_cast<double>(row.absolute_scale_sum)
            << ", \"effective_scales\": "
            << static_cast<double>(row.effective_scales)
            << ", \"signed_scale_alignment\": "
            << static_cast<double>(row.signed_scale_alignment)
            << ", \"absolute_fraction\": "
            << static_cast<double>(row.absolute_fraction) << '}'
            << (index + 1 == count ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"exact_reconstruction\": "
        << (report.exact_reconstruction ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_projective_sum_proved\": false,\n"
        << "  \"remaining_requirement\": \"derive a summable bound across primitive signature shapes after each fixed projective ray is closed by plane-sphere incidence\",\n"
        << "  \"finite_ledger_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldRemainderProjectiveReport
LocalSldRemainderProjectiveLedger::analyze(
    const LocalSldRemainderSignatureReport& signatures) {
    LocalSldRemainderProjectiveReport report;
    report.cutoff = signatures.cutoff;
    report.signature_count = signatures.signature_count;
    report.excludes_signature_123 = signatures.excludes_signature_123;
    report.excludes_triple_family = signatures.excludes_triple_family;
    report.expected_power_one_total = signatures.power_one_total;
    std::map<Shape, Accumulator> groups;
    for (const LocalSldRemainderSignatureEntry& signature :
         signatures.signatures) {
        const auto [shape, scale] = primitive_shape(
            signature.squared_lengths);
        LocalSldRemainderProjectiveEntry& entry = groups[shape].entry;
        entry.primitive_squared_lengths = shape;
        ++entry.scale_count;
        entry.interactions += signature.interactions;
        entry.minimum_scale = entry.minimum_scale == 0
            ? scale : std::min(entry.minimum_scale, scale);
        entry.maximum_scale = std::max(entry.maximum_scale, scale);
        entry.bracket_total += signature.total;
        entry.power_one_total += signature.power_one_ratio;
        entry.absolute_scale_sum += std::abs(signature.power_one_ratio);
        entry.squared_scale_sum +=
            signature.power_one_ratio * signature.power_one_ratio;
    }
    report.shapes.reserve(groups.size());
    for (auto& [shape, accumulator] : groups) {
        static_cast<void>(shape);
        LocalSldRemainderProjectiveEntry& entry = accumulator.entry;
        if (entry.squared_scale_sum > 0.0L) {
            entry.effective_scales = entry.absolute_scale_sum *
                entry.absolute_scale_sum / entry.squared_scale_sum;
        }
        if (entry.absolute_scale_sum > 0.0L) {
            entry.signed_scale_alignment =
                std::abs(entry.power_one_total) /
                entry.absolute_scale_sum;
        }
        report.reconstructed_power_one_total += entry.power_one_total;
        report.absolute_projective_sum +=
            std::abs(entry.power_one_total);
        report.squared_projective_sum +=
            entry.power_one_total * entry.power_one_total;
        report.shapes.push_back(entry);
    }
    report.projective_shape_count = report.shapes.size();
    if (report.squared_projective_sum > 0.0L) {
        report.effective_projective_shapes =
            report.absolute_projective_sum *
            report.absolute_projective_sum /
            report.squared_projective_sum;
    }
    if (report.absolute_projective_sum > 0.0L) {
        report.signed_projective_alignment =
            std::abs(report.reconstructed_power_one_total) /
            report.absolute_projective_sum;
        for (LocalSldRemainderProjectiveEntry& entry : report.shapes) {
            entry.absolute_fraction = std::abs(entry.power_one_total) /
                report.absolute_projective_sum;
            report.dominant_projective_fraction = std::max(
                report.dominant_projective_fraction,
                entry.absolute_fraction);
        }
    }
    std::sort(
        report.shapes.begin(), report.shapes.end(),
        [](const auto& left, const auto& right) {
            return std::abs(left.power_one_total) >
                std::abs(right.power_one_total);
        });
    report.reconstruction_error = relative_error(
        report.reconstructed_power_one_total,
        report.expected_power_one_total);
    report.exact_reconstruction = signatures.exact_reconstruction &&
        report.reconstruction_error < 1e-13L;
    return report;
}

LocalSldRemainderProjectiveCliOptions
LocalSldRemainderProjectiveCli::parse(
    int argc, char** argv, int first) {
    LocalSldRemainderProjectiveCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--top") {
            options.top = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--exclude-123") {
            options.exclude_signature_123 = true;
        } else if (name == "--exclude-triple-family") {
            options.exclude_triple_family = true;
        } else {
            throw std::invalid_argument(
                "unknown remainder-projective option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.top < 1 || options.threads < 1) {
        throw std::invalid_argument(
            "remainder-projective requires state, certificate, positive top and threads");
    }
    return options;
}

void LocalSldRemainderProjectiveCli::print_help(std::ostream& out) {
    out << "Local SLD remainder projective-signature options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON projective ledger\n"
        << "  --top N               store N strongest primitive shapes\n"
        << "  --threads N           parallel interaction workers\n"
        << "  --exclude-123         remove the fixed (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldRemainderProjectiveCli::run(
    const LocalSldRemainderProjectiveCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldRemainderSignatureReport signatures =
        LocalSldRemainderSignatureLedger::analyze(
            dynamics, state, options.threads,
            options.exclude_signature_123,
            options.exclude_triple_family);
    const LocalSldRemainderProjectiveReport report =
        LocalSldRemainderProjectiveLedger::analyze(signatures);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write remainder projective certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "remainder projective ledger cutoff=" << report.cutoff
        << " signatures=" << report.signature_count
        << " shapes=" << report.projective_shape_count
        << " effective="
        << static_cast<double>(report.effective_projective_shapes)
        << " dominant="
        << static_cast<double>(report.dominant_projective_fraction)
        << " power_one="
        << static_cast<double>(report.reconstructed_power_one_total)
        << " reconstruction="
        << static_cast<double>(report.reconstruction_error) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_reconstruction ? 0 : 2;
}

}  // namespace lemma
