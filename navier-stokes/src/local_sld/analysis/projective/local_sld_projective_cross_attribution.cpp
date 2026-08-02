#include "local_sld_projective_cross_attribution.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

void write_json(
    const LocalSldProjectiveCrossAttributionReport& report,
    const LocalSldProjectiveCrossAttributionCliOptions& options,
    std::ostream& output) {
    const std::size_t count = std::min(
        report.shapes.size(), static_cast<std::size_t>(options.top));
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-cross-attribution-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"projective_shape_count\": "
        << report.projective_shape_count << ",\n"
        << "  \"expected_cross_power_one\": "
        << static_cast<double>(report.expected_cross_power_one) << ",\n"
        << "  \"reconstructed_cross_power_one\": "
        << static_cast<double>(report.reconstructed_cross_power_one)
        << ",\n"
        << "  \"reconstruction_error\": "
        << static_cast<double>(report.reconstruction_error) << ",\n"
        << "  \"absolute_cross_attribution_sum\": "
        << static_cast<double>(report.absolute_cross_attribution_sum)
        << ",\n"
        << "  \"effective_cross_attribution_shapes\": "
        << static_cast<double>(report.effective_cross_attribution_shapes)
        << ",\n"
        << "  \"dominant_cross_attribution_fraction\": "
        << static_cast<double>(report.dominant_cross_attribution_fraction)
        << ",\n"
        << "  \"signed_cross_attribution_alignment\": "
        << static_cast<double>(report.signed_cross_attribution_alignment)
        << ",\n"
        << "  \"top_cross_attributions\": [\n";
    for (std::size_t index = 0; index < count; ++index) {
        const auto& row = report.shapes[index];
        output << "    {\"primitive_squared_lengths\": ["
            << row.primitive_squared_lengths[0] << ", "
            << row.primitive_squared_lengths[1] << ", "
            << row.primitive_squared_lengths[2]
            << "], \"full_attributed_power_one\": "
            << static_cast<double>(row.full_attributed_power_one)
            << ", \"self_power_one\": "
            << static_cast<double>(row.self_power_one)
            << ", \"cross_attributed_power_one\": "
            << static_cast<double>(row.cross_attributed_power_one)
            << ", \"absolute_fraction\": "
            << static_cast<double>(row.absolute_fraction) << '}'
            << (index + 1 == count ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"exact_reconstruction\": "
        << (report.exact_reconstruction ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_cross_bound_proved\": false,\n"
        << "  \"finite_ledger_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveCrossAttributionReport
LocalSldProjectiveCrossAttribution::analyze(
    const LocalSldRemainderProjectiveReport& full_attribution,
    const LocalSldProjectiveQuarticCrossReport& self_cross) {
    if (full_attribution.cutoff != self_cross.cutoff ||
        full_attribution.excludes_signature_123 !=
            self_cross.excludes_signature_123 ||
        full_attribution.excludes_triple_family !=
            self_cross.excludes_triple_family) {
        throw std::invalid_argument(
            "projective cross attribution report mismatch");
    }
    std::map<Shape, SpectralReal> full_by_shape;
    for (const auto& entry : full_attribution.shapes) {
        full_by_shape[entry.primitive_squared_lengths] =
            entry.power_one_total;
    }
    std::map<Shape, SpectralReal> self_by_shape;
    for (const auto& entry : self_cross.shapes) {
        self_by_shape[entry.primitive_squared_lengths] =
            entry.power_one_self;
    }
    LocalSldProjectiveCrossAttributionReport report;
    report.cutoff = self_cross.cutoff;
    report.expected_cross_power_one = self_cross.cross_power_one;
    for (const auto& [shape, full] : full_by_shape) {
        LocalSldProjectiveCrossAttributionEntry entry;
        entry.primitive_squared_lengths = shape;
        entry.full_attributed_power_one = full;
        entry.self_power_one = self_by_shape[shape];
        entry.cross_attributed_power_one =
            entry.full_attributed_power_one - entry.self_power_one;
        report.reconstructed_cross_power_one +=
            entry.cross_attributed_power_one;
        report.absolute_cross_attribution_sum +=
            std::abs(entry.cross_attributed_power_one);
        report.squared_cross_attribution_sum +=
            entry.cross_attributed_power_one *
            entry.cross_attributed_power_one;
        report.shapes.push_back(entry);
    }
    report.projective_shape_count = report.shapes.size();
    if (report.squared_cross_attribution_sum > 0.0L) {
        report.effective_cross_attribution_shapes =
            report.absolute_cross_attribution_sum *
            report.absolute_cross_attribution_sum /
            report.squared_cross_attribution_sum;
    }
    if (report.absolute_cross_attribution_sum > 0.0L) {
        report.signed_cross_attribution_alignment =
            std::abs(report.reconstructed_cross_power_one) /
            report.absolute_cross_attribution_sum;
        for (auto& entry : report.shapes) {
            entry.absolute_fraction =
                std::abs(entry.cross_attributed_power_one) /
                report.absolute_cross_attribution_sum;
            report.dominant_cross_attribution_fraction = std::max(
                report.dominant_cross_attribution_fraction,
                entry.absolute_fraction);
        }
    }
    std::sort(report.shapes.begin(), report.shapes.end(),
        [](const auto& left, const auto& right) {
            return std::abs(left.cross_attributed_power_one) >
                std::abs(right.cross_attributed_power_one);
        });
    report.reconstruction_error = std::abs(
        report.reconstructed_cross_power_one -
        report.expected_cross_power_one) /
        std::max({std::abs(report.expected_cross_power_one),
                  std::abs(full_attribution.expected_power_one_total),
                  std::abs(self_cross.diagonal_power_one),
                  1e-30L});
    report.exact_reconstruction =
        full_attribution.exact_reconstruction &&
        self_cross.exact_decomposition &&
        report.projective_shape_count == self_cross.projective_shape_count &&
        report.reconstruction_error < 1e-13L;
    return report;
}

LocalSldProjectiveCrossAttributionCliOptions
LocalSldProjectiveCrossAttributionCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveCrossAttributionCliOptions options;
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
                "unknown projective-cross-attribution option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.top < 1 || options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "projective-cross-attribution requires state, certificate, positive top, and threads");
    }
    return options;
}

void LocalSldProjectiveCrossAttributionCli::print_help(
    std::ostream& out) {
    out << "Projective cross-ray attribution options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON attribution\n"
        << "  --top N               reported primitive shapes\n"
        << "  --threads N           direct ledger workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveCrossAttributionCli::run(
    const LocalSldProjectiveCrossAttributionCliOptions& options,
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
    const LocalSldRemainderProjectiveReport full =
        LocalSldRemainderProjectiveLedger::analyze(signatures);
    const LocalSldProjectiveQuarticCrossReport self_cross =
        LocalSldProjectiveQuarticCrossLedger::analyze(
            dynamics, state, options.threads,
            options.exclude_signature_123,
            options.exclude_triple_family);
    const LocalSldProjectiveCrossAttributionReport report =
        LocalSldProjectiveCrossAttribution::analyze(full, self_cross);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective cross attribution certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective cross attribution cutoff=" << report.cutoff
        << " shapes=" << report.projective_shape_count
        << " effective="
        << static_cast<double>(report.effective_cross_attribution_shapes)
        << " dominant="
        << static_cast<double>(report.dominant_cross_attribution_fraction)
        << " alignment="
        << static_cast<double>(report.signed_cross_attribution_alignment)
        << " cross="
        << static_cast<double>(report.reconstructed_cross_power_one)
        << " error=" << static_cast<double>(report.reconstruction_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_reconstruction ? 0 : 2;
}

}  // namespace lemma
