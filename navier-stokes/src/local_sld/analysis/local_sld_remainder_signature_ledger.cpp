#include "local_sld_remainder_signature_ledger.hpp"

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
#include <vector>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {
namespace {

using Signature = std::array<SpectralInteger, 3>;

struct SelectedInteraction {
    InteractionIndex interaction{};
    std::size_t signature_index = 0;
};

struct Accumulator {
    std::size_t interactions = 0;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
    SpectralReal outer_pairing = 0.0L;
    SpectralReal advected_pairing = 0.0L;
    SpectralReal nested_pairing = 0.0L;
};

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "remainder signature Laplacian layout mismatch");
    }
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

void project_increment(
    SpectralIncrement& increment,
    const SpectralState& state) {
    if (increment.size() != state.waves.size()) {
        throw std::invalid_argument(
            "remainder signature projection layout mismatch");
    }
    for (std::size_t mode = 0; mode < increment.size(); ++mode) {
        increment[mode] = project_divergence_free(
            state.waves[mode], increment[mode]);
    }
}

SpectralReal relative_error(SpectralReal computed,
                            SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

Signature interaction_signature(
    const SpectralState& state,
    InteractionIndex interaction) {
    const auto [p, q, target] = interaction;
    Signature signature{
        norm_squared(state.waves[p]),
        norm_squared(state.waves[q]),
        norm_squared(state.waves[target])};
    std::sort(signature.begin(), signature.end());
    return signature;
}

SpectralReal raw_pairing(
    SpectralComplex coefficient,
    const ComplexVector& advected,
    const ComplexVector& cotangent) {
    SpectralReal result = 0.0L;
    for (std::size_t component = 0; component < 3; ++component) {
        result += std::real(
            std::conj(coefficient * advected[component]) *
            cotangent[component]);
    }
    return result;
}

void add(Accumulator& target, const Accumulator& source) {
    target.interactions += source.interactions;
    target.stretching += source.stretching;
    target.palinstrophy_cross += source.palinstrophy_cross;
    target.outer_pairing += source.outer_pairing;
    target.advected_pairing += source.advected_pairing;
    target.nested_pairing += source.nested_pairing;
}

void write_json(
    const LocalSldRemainderSignatureReport& report,
    const LocalSldRemainderSignatureCliOptions& options,
    std::ostream& output) {
    const std::size_t count = std::min(
        report.signatures.size(), static_cast<std::size_t>(options.top));
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-remainder-signature-ledger-v2\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"excludes_signature_123\": "
        << (report.excludes_signature_123 ? "true" : "false")
        << ",\n"
        << "  \"excludes_triple_family\": "
        << (report.excludes_triple_family ? "true" : "false")
        << ",\n"
        << "  \"selected_interactions\": "
        << report.selected_interactions << ",\n"
        << "  \"signature_count\": " << report.signature_count << ",\n"
        << "  \"reported_signature_count\": " << count << ",\n"
        << "  \"enstrophy\": "
        << static_cast<double>(report.enstrophy) << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(report.palinstrophy) << ",\n"
        << "  \"target_scale\": "
        << static_cast<double>(report.target_scale) << ",\n"
        << "  \"full_bracket\": "
        << static_cast<double>(report.full_bracket) << ",\n"
        << "  \"full_target_ratio\": "
        << static_cast<double>(report.full_target_ratio) << ",\n"
        << "  \"common_full_stretching\": "
        << static_cast<double>(report.common_full_stretching) << ",\n"
        << "  \"power_one_scale\": "
        << static_cast<double>(report.power_one_scale) << ",\n"
        << "  \"power_one_total\": "
        << static_cast<double>(report.power_one_total) << ",\n"
        << "  \"power_one_reconstruction_error\": "
        << static_cast<double>(report.power_one_reconstruction_error)
        << ",\n"
        << "  \"reconstructed_stretching\": "
        << static_cast<double>(report.reconstructed_stretching) << ",\n"
        << "  \"reconstructed_palinstrophy_cross\": "
        << static_cast<double>(
               report.reconstructed_palinstrophy_cross) << ",\n"
        << "  \"reconstructed_terms\": {\n"
        << "    \"outer_square\": "
        << static_cast<double>(report.reconstructed_outer_square) << ",\n"
        << "    \"advected_commutator\": "
        << static_cast<double>(
               report.reconstructed_advected_commutator) << ",\n"
        << "    \"enstrophy_normalization\": "
        << static_cast<double>(
               report.reconstructed_enstrophy_normalization) << ",\n"
        << "    \"palinstrophy_normalization\": "
        << static_cast<double>(
               report.reconstructed_palinstrophy_normalization) << ",\n"
        << "    \"advecting_nested\": "
        << static_cast<double>(
               report.reconstructed_advecting_nested) << ",\n"
        << "    \"total\": "
        << static_cast<double>(report.reconstructed_bracket) << "\n"
        << "  },\n"
        << "  \"stretching_reconstruction_error\": "
        << static_cast<double>(
               report.stretching_reconstruction_error) << ",\n"
        << "  \"palinstrophy_cross_reconstruction_error\": "
        << static_cast<double>(
               report.palinstrophy_cross_reconstruction_error) << ",\n"
        << "  \"bracket_reconstruction_error\": "
        << static_cast<double>(report.bracket_reconstruction_error)
        << ",\n"
        << "  \"absolute_contribution_sum\": "
        << static_cast<double>(report.absolute_contribution_sum) << ",\n"
        << "  \"squared_contribution_sum\": "
        << static_cast<double>(report.squared_contribution_sum) << ",\n"
        << "  \"effective_contributing_signatures\": "
        << static_cast<double>(
               report.effective_contributing_signatures) << ",\n"
        << "  \"dominant_absolute_fraction\": "
        << static_cast<double>(report.dominant_absolute_fraction) << ",\n"
        << "  \"signed_cancellation_ratio\": "
        << static_cast<double>(report.signed_cancellation_ratio) << ",\n"
        << "  \"signed_amplification\": "
        << static_cast<double>(report.signed_amplification) << ",\n"
        << "  \"top_signatures\": [\n";
    for (std::size_t index = 0; index < count; ++index) {
        const LocalSldRemainderSignatureEntry& row =
            report.signatures[index];
        output << "    {\"squared_lengths\": ["
            << row.squared_lengths[0] << ", "
            << row.squared_lengths[1] << ", "
            << row.squared_lengths[2] << "], \"interactions\": "
            << row.interactions
            << ", \"stretching\": "
            << static_cast<double>(row.stretching)
            << ", \"palinstrophy_cross\": "
            << static_cast<double>(row.palinstrophy_cross)
            << ", \"outer_square\": "
            << static_cast<double>(row.outer_square)
            << ", \"advected_commutator\": "
            << static_cast<double>(row.advected_commutator)
            << ", \"enstrophy_normalization\": "
            << static_cast<double>(row.enstrophy_normalization)
            << ", \"palinstrophy_normalization\": "
            << static_cast<double>(row.palinstrophy_normalization)
            << ", \"advecting_nested\": "
            << static_cast<double>(row.advecting_nested)
            << ", \"total\": " << static_cast<double>(row.total)
            << ", \"absolute_fraction\": "
            << static_cast<double>(row.absolute_fraction)
            << ", \"target_ratio\": "
            << static_cast<double>(row.target_ratio)
            << ", \"power_one_ratio\": "
            << static_cast<double>(row.power_one_ratio) << '}'
            << (index + 1 == count ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"exact_reconstruction\": "
        << (report.exact_reconstruction ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_bound_proved\": false,\n"
        << "  \"finite_ledger_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldRemainderSignatureReport
LocalSldRemainderSignatureLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (threads < 1) {
        throw std::invalid_argument(
            "remainder signature threads must be positive");
    }
    const TriadSelection selection = exclude_triple_family
        ? (exclude_signature_123
               ? TriadSelection::
                     local_without_equal_low_double_triple_and_signature(
                         1, 2, 3)
               : TriadSelection::local_without_equal_low_double_triple())
        : (exclude_signature_123
               ? TriadSelection::
                     local_without_equal_low_doubling_and_signature(
                         1, 2, 3)
               : TriadSelection::local_without_equal_low_doubling());
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(dynamics, selection).evaluate(state);
    const LocalQuarticClosureObjectiveValue common =
        LocalQuarticClosureObjective(
            dynamics, TriadPartition::local).evaluate(state);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const SpectralIncrement aau = laplacian_weight(state, au);
    const SpectralIncrement b =
        dynamics.advection_direct_partition(state, selection);
    const SpectralIncrement ab = laplacian_weight(state, b);
    const SpectralIncrement transported_au =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, au, selection);
    const BilinearAdvectionCotangents nested_cotangents =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, state.velocity, state.velocity, au, selection);

    SpectralIncrement projected_au = au;
    SpectralIncrement projected_aau = aau;
    SpectralIncrement projected_ab = ab;
    SpectralIncrement projected_transported_au = transported_au;
    SpectralIncrement projected_nested = nested_cotangents.advecting;
    project_increment(projected_au, state);
    project_increment(projected_aau, state);
    project_increment(projected_ab, state);
    project_increment(projected_transported_au, state);
    project_increment(projected_nested, state);

    std::map<Signature, std::size_t> signature_indices;
    std::vector<LocalSldRemainderSignatureEntry> entries;
    std::vector<SelectedInteraction> selected;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (!TriadPartitioner::includes(state, interaction, selection)) {
            continue;
        }
        const Signature signature = interaction_signature(state, interaction);
        auto [position, inserted] = signature_indices.try_emplace(
            signature, entries.size());
        if (inserted) {
            LocalSldRemainderSignatureEntry entry;
            entry.squared_lengths = signature;
            entries.push_back(entry);
        }
        selected.push_back({interaction, position->second});
    }

    int worker_count = 1;
#ifdef NS_HAVE_OPENMP
    worker_count = std::max(
        1, std::min(threads, omp_get_max_threads()));
#endif
    std::vector<std::vector<Accumulator>> partials(
        static_cast<std::size_t>(worker_count),
        std::vector<Accumulator>(entries.size()));
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count) \
    if(worker_count > 1)
#endif
    for (std::ptrdiff_t selected_index = 0;
         selected_index < static_cast<std::ptrdiff_t>(selected.size());
         ++selected_index) {
        int worker = 0;
#ifdef NS_HAVE_OPENMP
        worker = omp_get_thread_num();
#endif
        const SelectedInteraction& item = selected[
            static_cast<std::size_t>(selected_index)];
        const auto [p, q, target] = item.interaction;
        const SpectralComplex coefficient = imaginary_unit * wave_dot(
            state.waves[q], state.velocity[p]);
        const ComplexVector& advected = state.velocity[q];
        Accumulator& value = partials[
            static_cast<std::size_t>(worker)][item.signature_index];
        ++value.interactions;
        value.stretching += raw_pairing(
            coefficient, advected, projected_au[target]);
        value.palinstrophy_cross += raw_pairing(
            coefficient, advected, projected_aau[target]);
        value.outer_pairing += raw_pairing(
            coefficient, advected, projected_ab[target]);
        value.advected_pairing += raw_pairing(
            coefficient, advected, projected_transported_au[target]);
        value.nested_pairing += raw_pairing(
            coefficient, advected, projected_nested[target]);
    }

    std::vector<Accumulator> totals(entries.size());
    for (const std::vector<Accumulator>& partial : partials) {
        for (std::size_t index = 0; index < entries.size(); ++index) {
            add(totals[index], partial[index]);
        }
    }

    LocalSldRemainderSignatureReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.selected_interactions = selected.size();
    report.signature_count = entries.size();
    report.threads = worker_count;
    report.excludes_signature_123 = exclude_signature_123;
    report.excludes_triple_family = exclude_triple_family;
    report.enstrophy = full.enstrophy;
    report.palinstrophy = full.palinstrophy;
    report.target_scale = full.lqc3_target_scale;
    report.full_bracket = full.signed_two_entry_bracket;
    report.full_target_ratio = full.lqc3_target_ratio;
    report.common_full_stretching = common.signed_stretching;
    report.power_one_scale = common.signed_stretching /
        (common.enstrophy * common.enstrophy *
         common.palinstrophy * common.palinstrophy);
    for (std::size_t index = 0; index < entries.size(); ++index) {
        LocalSldRemainderSignatureEntry& entry = entries[index];
        const Accumulator& value = totals[index];
        entry.interactions = value.interactions;
        entry.stretching = value.stretching;
        entry.palinstrophy_cross = value.palinstrophy_cross;
        entry.outer_square = -value.outer_pairing;
        entry.advected_commutator = value.advected_pairing;
        entry.enstrophy_normalization = entry.stretching *
            full.signed_stretching / (2.0L * full.enstrophy);
        entry.palinstrophy_normalization = 3.0L * entry.stretching *
            full.palinstrophy_cross / (2.0L * full.palinstrophy);
        entry.advecting_nested = -value.nested_pairing;
        entry.total = entry.outer_square + entry.advected_commutator +
            entry.enstrophy_normalization +
            entry.palinstrophy_normalization +
            entry.advecting_nested;
        if (report.target_scale > 0.0L) {
            entry.target_ratio = entry.total / report.target_scale;
        }
        entry.power_one_ratio = entry.total * report.power_one_scale;
        report.reconstructed_stretching += entry.stretching;
        report.reconstructed_palinstrophy_cross +=
            entry.palinstrophy_cross;
        report.reconstructed_outer_square += entry.outer_square;
        report.reconstructed_advected_commutator +=
            entry.advected_commutator;
        report.reconstructed_enstrophy_normalization +=
            entry.enstrophy_normalization;
        report.reconstructed_palinstrophy_normalization +=
            entry.palinstrophy_normalization;
        report.reconstructed_advecting_nested +=
            entry.advecting_nested;
        report.reconstructed_bracket += entry.total;
        report.power_one_total += entry.power_one_ratio;
        report.absolute_contribution_sum += std::abs(entry.total);
        report.squared_contribution_sum += entry.total * entry.total;
    }
    if (report.absolute_contribution_sum > 0.0L) {
        for (LocalSldRemainderSignatureEntry& entry : entries) {
            entry.absolute_fraction = std::abs(entry.total) /
                report.absolute_contribution_sum;
            report.dominant_absolute_fraction = std::max(
                report.dominant_absolute_fraction,
                entry.absolute_fraction);
        }
        report.signed_cancellation_ratio = std::abs(
            report.reconstructed_bracket) /
            report.absolute_contribution_sum;
    }
    if (report.squared_contribution_sum > 0.0L) {
        report.effective_contributing_signatures =
            report.absolute_contribution_sum *
            report.absolute_contribution_sum /
            report.squared_contribution_sum;
        report.signed_amplification = std::abs(
            report.reconstructed_bracket) /
            std::sqrt(report.squared_contribution_sum);
    }
    report.stretching_reconstruction_error = relative_error(
        report.reconstructed_stretching, full.signed_stretching);
    report.palinstrophy_cross_reconstruction_error = relative_error(
        report.reconstructed_palinstrophy_cross,
        full.palinstrophy_cross);
    report.bracket_reconstruction_error = relative_error(
        report.reconstructed_bracket, report.full_bracket);
    report.power_one_reconstruction_error = relative_error(
        report.power_one_total,
        report.full_bracket * report.power_one_scale);
    report.exact_reconstruction =
        report.stretching_reconstruction_error < 1e-13L &&
        report.palinstrophy_cross_reconstruction_error < 1e-13L &&
        report.bracket_reconstruction_error < 1e-13L;
    report.exact_reconstruction = report.exact_reconstruction &&
        report.power_one_reconstruction_error < 1e-13L;
    std::sort(entries.begin(), entries.end(),
              [](const LocalSldRemainderSignatureEntry& left,
                 const LocalSldRemainderSignatureEntry& right) {
                  return std::abs(left.total) > std::abs(right.total);
              });
    report.signatures = std::move(entries);
    return report;
}

LocalSldRemainderSignatureCliOptions
LocalSldRemainderSignatureCli::parse(
    int argc, char** argv, int first) {
    LocalSldRemainderSignatureCliOptions options;
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
                "unknown remainder-signature option: " + name);
        }
    }
    if (options.state_path.empty() ||
        options.certificate_path.empty() ||
        options.top < 1 || options.threads < 1) {
        throw std::invalid_argument(
            "remainder-signature-ledger requires state, certificate, positive top and threads");
    }
    return options;
}

void LocalSldRemainderSignatureCli::print_help(std::ostream& out) {
    out << "Local SLD remainder signature-ledger options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON ledger\n"
        << "  --top N               store N strongest signature rows\n"
        << "  --exclude-123         also remove the fixed (1,2,3) signature\n"
        << "  --exclude-triple-family  also remove every (m,m,3m) signature\n"
        << "  --threads N           parallel interaction workers\n";
}

int LocalSldRemainderSignatureCli::run(
    const LocalSldRemainderSignatureCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldRemainderSignatureReport report =
        LocalSldRemainderSignatureLedger::analyze(
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
            "cannot write remainder signature ledger");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "remainder signature ledger cutoff=" << report.cutoff
        << " signatures=" << report.signature_count
        << " interactions=" << report.selected_interactions
        << " workers=" << report.threads
        << " lqc3=" << static_cast<double>(report.full_target_ratio)
        << " effective="
        << static_cast<double>(
               report.effective_contributing_signatures)
        << " dominant_fraction="
        << static_cast<double>(report.dominant_absolute_fraction)
        << " cancellation="
        << static_cast<double>(report.signed_cancellation_ratio)
        << " reconstruction="
        << static_cast<double>(report.bracket_reconstruction_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_reconstruction ? 0 : 2;
}

}  // namespace lemma
