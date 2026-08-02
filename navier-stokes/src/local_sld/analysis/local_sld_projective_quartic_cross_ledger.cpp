#include "local_sld_projective_quartic_cross_ledger.hpp"

#include "projective_advection_decomposition.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {
namespace {

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    SpectralIncrement result = source;
    if (result.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective quartic Laplacian layout mismatch");
    }
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective quartic pairing layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal relative_error(SpectralReal left, SpectralReal right) {
    return std::abs(left - right) /
        std::max({std::abs(left), std::abs(right), 1e-30L});
}

TriadSelection selection_for(
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (exclude_triple_family) {
        return exclude_signature_123
            ? TriadSelection::
                  local_without_equal_low_double_triple_and_signature(
                      1, 2, 3)
            : TriadSelection::local_without_equal_low_double_triple();
    }
    return exclude_signature_123
        ? TriadSelection::
              local_without_equal_low_doubling_and_signature(1, 2, 3)
        : TriadSelection::local_without_equal_low_doubling();
}

void write_json(
    const LocalSldProjectiveQuarticCrossReport& report,
    const LocalSldProjectiveQuarticCrossCliOptions& options,
    std::ostream& output) {
    const std::size_t count = std::min(
        report.shapes.size(), static_cast<std::size_t>(options.top));
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-quartic-cross-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"projective_shape_count\": "
        << report.projective_shape_count << ",\n"
        << "  \"excludes_signature_123\": "
        << (report.excludes_signature_123 ? "true" : "false")
        << ",\n"
        << "  \"excludes_triple_family\": "
        << (report.excludes_triple_family ? "true" : "false")
        << ",\n"
        << "  \"enstrophy\": " << static_cast<double>(report.enstrophy)
        << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(report.palinstrophy) << ",\n"
        << "  \"full_stretching\": "
        << static_cast<double>(report.full_stretching) << ",\n"
        << "  \"full_bracket\": "
        << static_cast<double>(report.full_bracket) << ",\n"
        << "  \"diagonal_bracket\": "
        << static_cast<double>(report.diagonal_bracket) << ",\n"
        << "  \"cross_bracket\": "
        << static_cast<double>(report.cross_bracket) << ",\n"
        << "  \"components\": {\n"
        << "    \"full_outer_commutator\": "
        << static_cast<double>(report.full_outer_commutator) << ",\n"
        << "    \"full_advecting_nested\": "
        << static_cast<double>(report.full_advecting_nested) << ",\n"
        << "    \"full_enstrophy_normalization\": "
        << static_cast<double>(report.full_enstrophy_normalization)
        << ",\n"
        << "    \"full_palinstrophy_normalization\": "
        << static_cast<double>(report.full_palinstrophy_normalization)
        << ",\n"
        << "    \"diagonal_outer_commutator\": "
        << static_cast<double>(report.diagonal_outer_commutator) << ",\n"
        << "    \"diagonal_advecting_nested\": "
        << static_cast<double>(report.diagonal_advecting_nested) << ",\n"
        << "    \"diagonal_enstrophy_normalization\": "
        << static_cast<double>(report.diagonal_enstrophy_normalization)
        << ",\n"
        << "    \"diagonal_palinstrophy_normalization\": "
        << static_cast<double>(report.diagonal_palinstrophy_normalization)
        << ",\n"
        << "    \"cross_outer_commutator\": "
        << static_cast<double>(report.cross_outer_commutator) << ",\n"
        << "    \"cross_advecting_nested\": "
        << static_cast<double>(report.cross_advecting_nested) << ",\n"
        << "    \"cross_enstrophy_normalization\": "
        << static_cast<double>(report.cross_enstrophy_normalization)
        << ",\n"
        << "    \"cross_palinstrophy_normalization\": "
        << static_cast<double>(report.cross_palinstrophy_normalization)
        << "\n  },\n"
        << "  \"power_one_scale\": "
        << static_cast<double>(report.power_one_scale) << ",\n"
        << "  \"full_power_one\": "
        << static_cast<double>(report.full_power_one) << ",\n"
        << "  \"diagonal_power_one\": "
        << static_cast<double>(report.diagonal_power_one) << ",\n"
        << "  \"cross_power_one\": "
        << static_cast<double>(report.cross_power_one) << ",\n"
        << "  \"cross_component_power_one\": {\n"
        << "    \"outer_commutator\": "
        << static_cast<double>(
               report.cross_outer_commutator_power_one) << ",\n"
        << "    \"advecting_nested\": "
        << static_cast<double>(
               report.cross_advecting_nested_power_one) << ",\n"
        << "    \"enstrophy_normalization\": "
        << static_cast<double>(
               report.cross_enstrophy_normalization_power_one) << ",\n"
        << "    \"palinstrophy_normalization\": "
        << static_cast<double>(
               report.cross_palinstrophy_normalization_power_one)
        << "\n  },\n"
        << "  \"absolute_diagonal_power_one_sum\": "
        << static_cast<double>(report.absolute_diagonal_power_one_sum)
        << ",\n"
        << "  \"effective_diagonal_shapes\": "
        << static_cast<double>(report.effective_diagonal_shapes) << ",\n"
        << "  \"dominant_diagonal_fraction\": "
        << static_cast<double>(report.dominant_diagonal_fraction)
        << ",\n"
        << "  \"diagonal_signed_alignment\": "
        << static_cast<double>(report.diagonal_signed_alignment) << ",\n"
        << "  \"full_component_reconstruction_error\": "
        << static_cast<double>(report.full_component_reconstruction_error)
        << ",\n"
        << "  \"bracket_decomposition_error\": "
        << static_cast<double>(report.bracket_decomposition_error)
        << ",\n"
        << "  \"top_diagonal_shapes\": [\n";
    for (std::size_t index = 0; index < count; ++index) {
        const auto& row = report.shapes[index];
        output << "    {\"primitive_squared_lengths\": ["
            << row.primitive_squared_lengths[0] << ", "
            << row.primitive_squared_lengths[1] << ", "
            << row.primitive_squared_lengths[2]
            << "], \"interactions\": " << row.interactions
            << ", \"self_bracket\": "
            << static_cast<double>(row.self_bracket)
            << ", \"power_one_self\": "
            << static_cast<double>(row.power_one_self)
            << ", \"absolute_fraction\": "
            << static_cast<double>(row.absolute_fraction) << '}'
            << (index + 1 == count ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"exact_decomposition\": "
        << (report.exact_decomposition ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_cross_bound_proved\": false,\n"
        << "  \"finite_ledger_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveQuarticCrossReport
LocalSldProjectiveQuarticCrossLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective quartic cross threads must be 1..256");
    }
    const TriadSelection selection = selection_for(
        exclude_signature_123, exclude_triple_family);
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(dynamics, selection).evaluate(state);
    const LocalQuarticClosureObjectiveValue full_local =
        LocalQuarticClosureObjective(
            dynamics, TriadPartition::local).evaluate(state);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    std::vector<LocalSldProjectiveQuarticSelfEntry> entries(
        groups.size());
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) \
    if(threads > 1)
#endif
    for (std::ptrdiff_t group_index = 0;
         group_index < static_cast<std::ptrdiff_t>(groups.size());
         ++group_index) {
        const std::size_t index = static_cast<std::size_t>(group_index);
        const ProjectiveInteractionGroup& group = groups[index];
        LocalSldProjectiveQuarticSelfEntry& entry = entries[index];
        entry.primitive_squared_lengths =
            group.primitive_squared_lengths;
        entry.interactions = group.interactions.size();
        const SpectralIncrement b =
            ProjectiveAdvectionDecomposition::evaluate(state, group);
        const SpectralIncrement ab = laplacian_weight(state, b);
        const SpectralIncrement transported_au =
            ProjectiveAdvectionDecomposition::evaluate_bilinear(
                state, group, state.velocity, au);
        const SpectralIncrement b_advects_u =
            ProjectiveAdvectionDecomposition::evaluate_bilinear(
                state, group, b, state.velocity);
        const SpectralReal stretching = pairing(au, b);
        const SpectralReal cross = pairing(ab, au);
        entry.outer_square = -pairing(b, ab);
        entry.advected_commutator = pairing(b, transported_au);
        entry.advecting_nested = -pairing(au, b_advects_u);
        entry.enstrophy_normalization = stretching * stretching /
            (2.0L * selected.enstrophy);
        entry.palinstrophy_normalization =
            3.0L * stretching * cross /
            (2.0L * selected.palinstrophy);
        entry.self_bracket = entry.outer_square +
            entry.advected_commutator + entry.advecting_nested +
            entry.enstrophy_normalization +
            entry.palinstrophy_normalization;
    }

    LocalSldProjectiveQuarticCrossReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.threads = threads;
    report.excludes_signature_123 = exclude_signature_123;
    report.excludes_triple_family = exclude_triple_family;
    report.projective_shape_count = entries.size();
    report.enstrophy = selected.enstrophy;
    report.palinstrophy = selected.palinstrophy;
    report.full_stretching = full_local.signed_stretching;
    report.full_bracket = selected.signed_two_entry_bracket;
    report.full_outer_commutator =
        selected.negative_commutator_pairing;
    report.full_advecting_nested = selected.advecting_slot;
    report.full_enstrophy_normalization =
        selected.signed_stretching * selected.signed_stretching /
        (2.0L * selected.enstrophy);
    report.full_palinstrophy_normalization =
        3.0L * selected.signed_stretching *
        selected.palinstrophy_cross /
        (2.0L * selected.palinstrophy);
    if (report.enstrophy > 0.0L && report.palinstrophy > 0.0L) {
        report.power_one_scale = report.full_stretching /
            (report.enstrophy * report.enstrophy *
             report.palinstrophy * report.palinstrophy);
    }
    SpectralReal diagonal_square_sum = 0.0L;
    for (auto& entry : entries) {
        report.diagonal_outer_commutator +=
            entry.outer_square + entry.advected_commutator;
        report.diagonal_advecting_nested += entry.advecting_nested;
        report.diagonal_enstrophy_normalization +=
            entry.enstrophy_normalization;
        report.diagonal_palinstrophy_normalization +=
            entry.palinstrophy_normalization;
        report.diagonal_bracket += entry.self_bracket;
        entry.power_one_self =
            entry.self_bracket * report.power_one_scale;
        report.absolute_diagonal_power_one_sum +=
            std::abs(entry.power_one_self);
        diagonal_square_sum +=
            entry.power_one_self * entry.power_one_self;
    }
    report.cross_bracket = report.full_bracket - report.diagonal_bracket;
    report.cross_outer_commutator = report.full_outer_commutator -
        report.diagonal_outer_commutator;
    report.cross_advecting_nested = report.full_advecting_nested -
        report.diagonal_advecting_nested;
    report.cross_enstrophy_normalization =
        report.full_enstrophy_normalization -
        report.diagonal_enstrophy_normalization;
    report.cross_palinstrophy_normalization =
        report.full_palinstrophy_normalization -
        report.diagonal_palinstrophy_normalization;
    report.full_power_one = report.full_bracket * report.power_one_scale;
    report.diagonal_power_one =
        report.diagonal_bracket * report.power_one_scale;
    report.cross_power_one =
        report.cross_bracket * report.power_one_scale;
    report.cross_outer_commutator_power_one =
        report.cross_outer_commutator * report.power_one_scale;
    report.cross_advecting_nested_power_one =
        report.cross_advecting_nested * report.power_one_scale;
    report.cross_enstrophy_normalization_power_one =
        report.cross_enstrophy_normalization * report.power_one_scale;
    report.cross_palinstrophy_normalization_power_one =
        report.cross_palinstrophy_normalization * report.power_one_scale;
    if (diagonal_square_sum > 0.0L) {
        report.effective_diagonal_shapes =
            report.absolute_diagonal_power_one_sum *
            report.absolute_diagonal_power_one_sum /
            diagonal_square_sum;
    }
    if (report.absolute_diagonal_power_one_sum > 0.0L) {
        report.diagonal_signed_alignment =
            std::abs(report.diagonal_power_one) /
            report.absolute_diagonal_power_one_sum;
        for (auto& entry : entries) {
            entry.absolute_fraction = std::abs(entry.power_one_self) /
                report.absolute_diagonal_power_one_sum;
            report.dominant_diagonal_fraction = std::max(
                report.dominant_diagonal_fraction,
                entry.absolute_fraction);
        }
    }
    std::sort(entries.begin(), entries.end(),
        [](const auto& left, const auto& right) {
            return std::abs(left.power_one_self) >
                std::abs(right.power_one_self);
        });
    report.shapes = std::move(entries);
    const SpectralReal reconstructed_full =
        report.full_outer_commutator +
        report.full_advecting_nested +
        report.full_enstrophy_normalization +
        report.full_palinstrophy_normalization;
    report.full_component_reconstruction_error = relative_error(
        reconstructed_full, report.full_bracket);
    report.bracket_decomposition_error = relative_error(
        report.diagonal_bracket + report.cross_bracket,
        report.full_bracket);
    report.exact_decomposition =
        report.full_component_reconstruction_error < 1e-13L &&
        report.bracket_decomposition_error < 1e-13L;
    return report;
}

LocalSldProjectiveQuarticCrossCliOptions
LocalSldProjectiveQuarticCrossCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveQuarticCrossCliOptions options;
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
                "unknown projective-quartic-cross option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.top < 1 || options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "projective-quartic-cross requires state, certificate, positive top, and threads");
    }
    return options;
}

void LocalSldProjectiveQuarticCrossCli::print_help(std::ostream& out) {
    out << "Projective quartic self/cross decomposition options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON ledger\n"
        << "  --top N               reported diagonal projective shapes\n"
        << "  --threads N           parallel projective-shape workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveQuarticCrossCli::run(
    const LocalSldProjectiveQuarticCrossCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldProjectiveQuarticCrossReport report =
        LocalSldProjectiveQuarticCrossLedger::analyze(
        dynamics, state, options.threads,
        options.exclude_signature_123,
        options.exclude_triple_family);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective quartic cross certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective quartic cross cutoff=" << report.cutoff
        << " shapes=" << report.projective_shape_count
        << " full=" << static_cast<double>(report.full_power_one)
        << " diagonal="
        << static_cast<double>(report.diagonal_power_one)
        << " cross=" << static_cast<double>(report.cross_power_one)
        << " effective_diagonal="
        << static_cast<double>(report.effective_diagonal_shapes)
        << " error="
        << static_cast<double>(report.full_component_reconstruction_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_decomposition ? 0 : 2;
}

}  // namespace lemma
