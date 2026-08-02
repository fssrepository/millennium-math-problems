#include "local_sld_projective_core_tail_ledger.hpp"

#include "projective_advection_decomposition.hpp"
#include "projective_triad_geometry.hpp"
#include "projective_core_family.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

using Signature = LocalSldProjectiveCoreSignature;

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective core-tail Laplacian layout mismatch");
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

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective core-tail pairing layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal relative_error(
    SpectralReal computed,
    SpectralReal expected,
    SpectralReal reference_scale = 0.0L) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected),
                  std::abs(reference_scale), 1e-30L});
}

Signature parse_signature(const std::string& text) {
    Signature result{};
    std::stringstream stream(text);
    std::string token;
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (!std::getline(stream, token, ',')) {
            throw std::invalid_argument(
                "core signature must have three comma-separated integers");
        }
        result[index] = static_cast<SpectralInteger>(std::stoll(token));
    }
    if (std::getline(stream, token, ',')) {
        throw std::invalid_argument(
            "core signature must have exactly three entries");
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<Signature> canonical_core(
    std::vector<Signature> core) {
    if (core.empty()) {
        throw std::invalid_argument(
            "projective core-tail ledger requires a nonempty fixed core");
    }
    for (Signature& signature : core) {
        std::sort(signature.begin(), signature.end());
        const auto geometry = ProjectiveTriadGeometry::certify(
            1, signature);
        if (!geometry.triangle_feasible) {
            throw std::invalid_argument(
                "projective core contains an infeasible signature");
        }
    }
    std::sort(core.begin(), core.end());
    if (std::adjacent_find(core.begin(), core.end()) != core.end()) {
        throw std::invalid_argument(
            "projective core contains a duplicate signature");
    }
    return core;
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

void finalize_block(
    LocalSldProjectiveCoreTailBlock& block,
    SpectralReal power_one_scale) {
    block.bracket = block.outer_square +
        block.advected_commutator + block.advecting_nested +
        block.enstrophy_normalization +
        block.palinstrophy_normalization;
    block.power_one = block.bracket * power_one_scale;
}

void write_signature(std::ostream& out, const Signature& signature) {
    out << '[' << signature[0] << ", " << signature[1]
        << ", " << signature[2] << ']';
}

void write_block(
    std::ostream& out,
    const std::string& name,
    const LocalSldProjectiveCoreTailBlock& block,
    bool trailing_comma) {
    out << "  \"" << name << "\": {\n"
        << "    \"outer_square\": "
        << static_cast<double>(block.outer_square) << ",\n"
        << "    \"advected_commutator\": "
        << static_cast<double>(block.advected_commutator) << ",\n"
        << "    \"advecting_nested\": "
        << static_cast<double>(block.advecting_nested) << ",\n"
        << "    \"enstrophy_normalization\": "
        << static_cast<double>(block.enstrophy_normalization) << ",\n"
        << "    \"palinstrophy_normalization\": "
        << static_cast<double>(block.palinstrophy_normalization) << ",\n"
        << "    \"bracket\": "
        << static_cast<double>(block.bracket) << ",\n"
        << "    \"power_one\": "
        << static_cast<double>(block.power_one) << "\n"
        << "  }" << (trailing_comma ? ",\n" : "\n");
}

void write_json(
    const LocalSldProjectiveCoreTailReport& report,
    const LocalSldProjectiveCoreTailCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-core-tail-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"requested_core\": [";
    for (std::size_t index = 0; index < report.requested_core.size();
         ++index) {
        if (index != 0) {
            output << ", ";
        }
        write_signature(output, report.requested_core[index]);
    }
    output << "],\n"
        << "  \"active_core_shape_count\": "
        << report.active_core_shape_count << ",\n"
        << "  \"tail_shape_count\": " << report.tail_shape_count
        << ",\n"
        << "  \"core_interaction_count\": "
        << report.core_interaction_count << ",\n"
        << "  \"tail_interaction_count\": "
        << report.tail_interaction_count << ",\n"
        << "  \"enstrophy\": "
        << static_cast<double>(report.enstrophy) << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(report.palinstrophy) << ",\n"
        << "  \"full_local_stretching\": "
        << static_cast<double>(report.full_local_stretching) << ",\n"
        << "  \"selected_stretching\": "
        << static_cast<double>(report.selected_stretching) << ",\n"
        << "  \"core_stretching\": "
        << static_cast<double>(report.core_stretching) << ",\n"
        << "  \"tail_stretching\": "
        << static_cast<double>(report.tail_stretching) << ",\n"
        << "  \"power_one_scale\": "
        << static_cast<double>(report.power_one_scale) << ",\n"
        << "  \"selected_bracket\": "
        << static_cast<double>(report.selected_bracket) << ",\n"
        << "  \"selected_power_one\": "
        << static_cast<double>(report.selected_power_one) << ",\n";
    write_block(output, "core_internal", report.core, true);
    write_block(output, "core_tail", report.core_tail, true);
    write_block(output, "tail_internal", report.tail, true);
    output << "  \"stretching_partition_error\": "
        << static_cast<double>(report.stretching_partition_error)
        << ",\n"
        << "  \"palinstrophy_cross_partition_error\": "
        << static_cast<double>(
               report.palinstrophy_cross_partition_error) << ",\n"
        << "  \"bracket_partition_error\": "
        << static_cast<double>(report.bracket_partition_error) << ",\n"
        << "  \"exact_core_tail_decomposition\": "
        << (report.exact_core_tail_decomposition ? "true" : "false")
        << ",\n"
        << "  \"fixed_core_internal_bound_proved\": true,\n"
        << "  \"core_tail_bound_proved\": false,\n"
        << "  \"growing_tail_internal_bound_proved\": false,\n"
        << "  \"full_local_lemma_proved\": false,\n"
        << "  \"remaining_requirement\": \"bound the signed core-tail and growing-tail internal power-one sum uniformly as the cutoff grows\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveCoreTailReport
LocalSldProjectiveCoreTailLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    std::vector<LocalSldProjectiveCoreSignature> core,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective core-tail threads must be 1..256");
    }
    core = canonical_core(std::move(core));
    const std::set<Signature> core_set(core.begin(), core.end());
    const TriadSelection selection = selection_for(
        exclude_signature_123, exclude_triple_family);
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    std::vector<std::size_t> core_indices;
    std::vector<std::size_t> tail_indices;
    std::vector<Signature> active_core;
    for (std::size_t index = 0; index < groups.size(); ++index) {
        if (core_set.contains(groups[index].primitive_squared_lengths)) {
            core_indices.push_back(index);
            active_core.push_back(groups[index].primitive_squared_lengths);
        } else {
            tail_indices.push_back(index);
        }
    }
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(dynamics, selection).evaluate(state);
    const LocalQuarticClosureObjectiveValue full_local =
        LocalQuarticClosureObjective(
            dynamics, TriadPartition::local).evaluate(state);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    auto evaluate = [&](const std::vector<std::size_t>& indices,
                        const SpectralIncrement& advecting,
                        const SpectralIncrement& advected) {
        return ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, indices, advecting, advected, threads);
    };
    const SpectralIncrement b_core = evaluate(
        core_indices, state.velocity, state.velocity);
    const SpectralIncrement b_tail = evaluate(
        tail_indices, state.velocity, state.velocity);
    const SpectralIncrement ab_core = laplacian_weight(state, b_core);
    const SpectralIncrement ab_tail = laplacian_weight(state, b_tail);
    const SpectralIncrement c_core = evaluate(
        core_indices, state.velocity, au);
    const SpectralIncrement c_tail = evaluate(
        tail_indices, state.velocity, au);
    const SpectralIncrement nested_core_core = evaluate(
        core_indices, b_core, state.velocity);
    const SpectralIncrement nested_tail_tail = evaluate(
        tail_indices, b_tail, state.velocity);
    const SpectralIncrement nested_core_tail = evaluate(
        core_indices, b_tail, state.velocity);
    const SpectralIncrement nested_tail_core = evaluate(
        tail_indices, b_core, state.velocity);

    LocalSldProjectiveCoreTailReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.threads = threads;
    report.excludes_signature_123 = exclude_signature_123;
    report.excludes_triple_family = exclude_triple_family;
    report.requested_core = core;
    report.active_core = active_core;
    report.requested_core_shape_count = core.size();
    report.active_core_shape_count = core_indices.size();
    report.tail_shape_count = tail_indices.size();
    for (const std::size_t index : core_indices) {
        report.core_interaction_count += groups[index].interactions.size();
    }
    for (const std::size_t index : tail_indices) {
        report.tail_interaction_count += groups[index].interactions.size();
    }
    report.enstrophy = selected.enstrophy;
    report.palinstrophy = selected.palinstrophy;
    report.full_local_stretching = full_local.signed_stretching;
    report.selected_stretching = selected.signed_stretching;
    report.core_stretching = pairing(au, b_core);
    report.tail_stretching = pairing(au, b_tail);
    report.core_palinstrophy_cross = pairing(ab_core, au);
    report.tail_palinstrophy_cross = pairing(ab_tail, au);
    if (report.enstrophy > 0.0L && report.palinstrophy > 0.0L) {
        report.power_one_scale = report.full_local_stretching /
            (report.enstrophy * report.enstrophy *
             report.palinstrophy * report.palinstrophy);
    }
    report.selected_bracket = selected.signed_two_entry_bracket;
    report.selected_power_one =
        report.selected_bracket * report.power_one_scale;

    report.core.outer_square = -pairing(b_core, ab_core);
    report.core.advected_commutator = pairing(b_core, c_core);
    report.core.advecting_nested = -pairing(
        au, nested_core_core);
    report.core.enstrophy_normalization =
        report.core_stretching * report.core_stretching /
        (2.0L * report.enstrophy);
    report.core.palinstrophy_normalization =
        3.0L * report.core_stretching *
        report.core_palinstrophy_cross /
        (2.0L * report.palinstrophy);
    finalize_block(report.core, report.power_one_scale);

    report.tail.outer_square = -pairing(b_tail, ab_tail);
    report.tail.advected_commutator = pairing(b_tail, c_tail);
    report.tail.advecting_nested = -pairing(
        au, nested_tail_tail);
    report.tail.enstrophy_normalization =
        report.tail_stretching * report.tail_stretching /
        (2.0L * report.enstrophy);
    report.tail.palinstrophy_normalization =
        3.0L * report.tail_stretching *
        report.tail_palinstrophy_cross /
        (2.0L * report.palinstrophy);
    finalize_block(report.tail, report.power_one_scale);

    report.core_tail.outer_square =
        -pairing(b_core, ab_tail) - pairing(b_tail, ab_core);
    report.core_tail.advected_commutator =
        pairing(b_core, c_tail) + pairing(b_tail, c_core);
    report.core_tail.advecting_nested = -pairing(
        au, nested_core_tail) - pairing(au, nested_tail_core);
    report.core_tail.enstrophy_normalization =
        report.core_stretching * report.tail_stretching /
        report.enstrophy;
    report.core_tail.palinstrophy_normalization =
        3.0L *
        (report.core_stretching *
             report.tail_palinstrophy_cross +
         report.tail_stretching *
             report.core_palinstrophy_cross) /
        (2.0L * report.palinstrophy);
    finalize_block(report.core_tail, report.power_one_scale);

    report.stretching_partition_error = relative_error(
        report.core_stretching + report.tail_stretching,
        report.selected_stretching,
        report.full_local_stretching);
    report.palinstrophy_cross_partition_error = relative_error(
        report.core_palinstrophy_cross +
            report.tail_palinstrophy_cross,
        selected.palinstrophy_cross);
    report.bracket_partition_error = relative_error(
        report.core.bracket + report.core_tail.bracket +
            report.tail.bracket,
        report.selected_bracket);
    report.exact_core_tail_decomposition =
        report.stretching_partition_error < 1e-13L &&
        report.palinstrophy_cross_partition_error < 1e-13L &&
        report.bracket_partition_error < 1e-13L;
    report.fixed_core_internal_bound_proved = true;
    return report;
}

LocalSldProjectiveCoreTailCliOptions
LocalSldProjectiveCoreTailCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveCoreTailCliOptions options;
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
        } else if (name == "--core-signature") {
            options.core.push_back(parse_signature(next(index, name)));
        } else if (name == "--core-max-height") {
            options.core_maximum_height =
                static_cast<SpectralInteger>(std::stoll(next(index, name)));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--exclude-123") {
            options.exclude_signature_123 = true;
        } else if (name == "--exclude-triple-family") {
            options.exclude_triple_family = true;
        } else {
            throw std::invalid_argument(
                "unknown projective-core-tail option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        (options.core.empty() && options.core_maximum_height == 0) ||
        (!options.core.empty() && options.core_maximum_height != 0) ||
        options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "projective-core-tail requires state, certificate, exactly one core definition, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveCoreTailCli::print_help(std::ostream& out) {
    out << "Projective fixed-core/growing-tail quartet options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON ledger\n"
        << "  --core-signature A,B,C  add one fixed primitive core ray; repeatable\n"
        << "  --core-max-height H   use every primitive feasible ray with max(a,b,c)<=H\n"
        << "  --threads N           parallel projective interaction workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveCoreTailCli::run(
    const LocalSldProjectiveCoreTailCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const std::vector<LocalSldProjectiveCoreSignature> core =
        options.core_maximum_height > 0
        ? ProjectiveCoreFamily::through_maximum_height(
              options.core_maximum_height)
        : options.core;
    const LocalSldProjectiveCoreTailReport report =
        LocalSldProjectiveCoreTailLedger::analyze(
            dynamics, state, core, options.threads,
            options.exclude_signature_123,
            options.exclude_triple_family);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective core-tail certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective core-tail cutoff=" << report.cutoff
        << " core_shapes=" << report.active_core_shape_count
        << " tail_shapes=" << report.tail_shape_count
        << " power_one(core,mixed,tail)=("
        << static_cast<double>(report.core.power_one) << ','
        << static_cast<double>(report.core_tail.power_one) << ','
        << static_cast<double>(report.tail.power_one) << ") error="
        << static_cast<double>(report.bracket_partition_error) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_core_tail_decomposition ? 0 : 2;
}

}  // namespace lemma
