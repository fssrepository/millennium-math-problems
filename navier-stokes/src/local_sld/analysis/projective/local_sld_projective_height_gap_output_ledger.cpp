#include "local_sld_projective_height_gap_output_ledger.hpp"

#include "local_sld_triad_selection.hpp"
#include "projective_advection_decomposition.hpp"
#include "projective_height_shell_partition.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace lemma {
namespace {

struct ShellField {
    std::size_t shape_count = 0;
    std::size_t interaction_count = 0;
    SpectralIncrement ab;
};

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "height-gap output Laplacian layout mismatch");
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

ShellField evaluate_shell(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    const ProjectiveHeightShellGroup& shell,
    int threads) {
    ShellField result;
    result.shape_count = shell.group_indices.size();
    result.interaction_count = shell.interaction_count;
    result.ab = laplacian_weight(
        state,
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, shell.group_indices,
            state.velocity, state.velocity, threads));
    return result;
}

SpectralReal norm2(const ComplexVector& value) {
    return std::real(dot_hermitian(value, value));
}

SpectralReal relative_error(
    SpectralReal computed,
    SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

void append_mode(
    LocalSldProjectiveHeightGapOutputMode& target,
    const ComplexVector& first,
    const ComplexVector& second) {
    const SpectralComplex value = dot_hermitian(first, second);
    target.first_h2_norm2 += norm2(first);
    target.second_h2_norm2 += norm2(second);
    target.pairing_real += std::real(value);
    target.pairing_imaginary += std::imag(value);
}

void normalize_mode(
    LocalSldProjectiveHeightGapOutputMode& mode,
    const LocalSldProjectiveHeightGapOutputReport& report) {
    mode.pairing_magnitude = std::abs(SpectralComplex(
        mode.pairing_real, mode.pairing_imaginary));
    const SpectralReal denominator = std::sqrt(std::max(
        mode.first_h2_norm2 * mode.second_h2_norm2, 0.0L));
    if (denominator > 0.0L) {
        mode.modal_correlation = mode.pairing_magnitude / denominator;
    }
    if (report.first_h2_norm2 > 0.0L) {
        mode.first_energy_fraction =
            mode.first_h2_norm2 / report.first_h2_norm2;
    }
    if (report.second_h2_norm2 > 0.0L) {
        mode.second_energy_fraction =
            mode.second_h2_norm2 / report.second_h2_norm2;
    }
    if (std::abs(report.gram_pairing) > 1e-30L) {
        mode.signed_gram_fraction =
            mode.pairing_real / report.gram_pairing;
    }
    if (report.absolute_modal_pairing > 0.0L) {
        mode.absolute_pairing_fraction =
            mode.pairing_magnitude / report.absolute_modal_pairing;
    }
}

}  // namespace

LocalSldProjectiveHeightGapOutputReport
LocalSldProjectiveHeightGapOutputLedger::analyze(
    const SpectralState& state,
    TriadSelection selection,
    int first_shell,
    int second_shell,
    std::size_t top_modes,
    int threads) {
    if (first_shell < 0 || second_shell <= first_shell) {
        throw std::invalid_argument(
            "height-gap output requires 0 <= first < second");
    }
    if (top_modes < 1 || top_modes > 1000 || threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "height-gap output requires top modes 1..1000 and threads 1..256");
    }
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    const auto partition = ProjectiveHeightShellPartition::build(groups);
    if (static_cast<std::size_t>(second_shell) >= partition.size() ||
        partition[static_cast<std::size_t>(first_shell)]
            .group_indices.empty() ||
        partition[static_cast<std::size_t>(second_shell)]
            .group_indices.empty()) {
        throw std::invalid_argument(
            "height-gap output shell is absent from the cutoff layout");
    }
    const ShellField first = evaluate_shell(
        state, groups, partition[static_cast<std::size_t>(first_shell)],
        threads);
    const ShellField second = evaluate_shell(
        state, groups, partition[static_cast<std::size_t>(second_shell)],
        threads);

    LocalSldProjectiveHeightGapOutputReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.first_shell = first_shell;
    report.second_shell = second_shell;
    report.shell_gap = second_shell - first_shell;
    report.first_shape_count = first.shape_count;
    report.second_shape_count = second.shape_count;
    report.first_interaction_count = first.interaction_count;
    report.second_interaction_count = second.interaction_count;

    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        report.first_h2_norm2 += norm2(first.ab[mode]);
        report.second_h2_norm2 += norm2(second.ab[mode]);
        report.gram_pairing += std::real(dot_hermitian(
            first.ab[mode], second.ab[mode]));
    }
    const SpectralReal denominator = std::sqrt(std::max(
        report.first_h2_norm2 * report.second_h2_norm2, 0.0L));
    if (denominator > 0.0L) {
        report.correlation = report.gram_pairing / denominator;
    }
    const SpectralReal first_threshold = std::max(
        report.first_h2_norm2 * 1e-24L, 1e-30L);
    const SpectralReal second_threshold = std::max(
        report.second_h2_norm2 * 1e-24L, 1e-30L);
    SpectralReal shared_first = 0.0L;
    SpectralReal shared_second = 0.0L;
    SpectralReal pairing_square_sum = 0.0L;
    std::vector<LocalSldProjectiveHeightGapOutputMode> modes;
    modes.reserve(state.waves.size() / 2U + 1U);
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        if (!is_positive_representative(wave)) {
            continue;
        }
        LocalSldProjectiveHeightGapOutputMode row;
        row.wave = wave;
        append_mode(row, first.ab[index], second.ab[index]);
        const auto negative = state.index.find(-wave);
        if (negative != state.index.end() && negative->second != index) {
            append_mode(
                row, first.ab[negative->second], second.ab[negative->second]);
        }
        row.pairing_magnitude = std::abs(SpectralComplex(
            row.pairing_real, row.pairing_imaginary));
        report.reconstructed_gram_pairing += row.pairing_real;
        report.absolute_modal_pairing += row.pairing_magnitude;
        pairing_square_sum += row.pairing_magnitude * row.pairing_magnitude;
        if (row.first_h2_norm2 > first_threshold &&
            row.second_h2_norm2 > second_threshold) {
            ++report.shared_conjugate_pair_count;
            shared_first += row.first_h2_norm2;
            shared_second += row.second_h2_norm2;
        }
        modes.push_back(std::move(row));
    }
    if (report.first_h2_norm2 > 0.0L) {
        report.first_shared_energy_fraction =
            shared_first / report.first_h2_norm2;
    }
    if (report.second_h2_norm2 > 0.0L) {
        report.second_shared_energy_fraction =
            shared_second / report.second_h2_norm2;
    }
    if (pairing_square_sum > 0.0L) {
        report.effective_shared_mode_count =
            report.absolute_modal_pairing * report.absolute_modal_pairing /
            pairing_square_sum;
    }
    for (auto& mode : modes) {
        normalize_mode(mode, report);
    }
    std::sort(
        modes.begin(), modes.end(),
        [](const auto& left, const auto& right) {
            return left.pairing_magnitude > right.pairing_magnitude;
        });
    const std::size_t retained = std::min(top_modes, modes.size());
    report.top_modes.assign(modes.begin(), modes.begin() + retained);
    for (const auto& mode : report.top_modes) {
        report.top_signed_gram_fraction += mode.signed_gram_fraction;
        report.top_absolute_pairing_fraction +=
            mode.absolute_pairing_fraction;
    }
    report.gram_reconstruction_error = relative_error(
        report.reconstructed_gram_pairing, report.gram_pairing);
    report.exact_gram_reconstruction =
        report.gram_reconstruction_error < 2e-12L;
    report.finite = report.exact_gram_reconstruction &&
        std::isfinite(report.correlation) &&
        std::isfinite(report.effective_shared_mode_count);
    return report;
}

LocalSldProjectiveHeightGapOutputCliOptions
LocalSldProjectiveHeightGapOutputCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightGapOutputCliOptions options;
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
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--first-shell") {
            options.first_shell = std::stoi(next(index, name));
        } else if (name == "--second-shell") {
            options.second_shell = std::stoi(next(index, name));
        } else if (name == "--top-modes") {
            options.top_modes = static_cast<std::size_t>(
                std::stoull(next(index, name)));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown height-gap output option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.first_shell < 0 ||
        options.second_shell <= options.first_shell ||
        options.second_shell > 62 || options.top_modes < 1 ||
        options.top_modes > 1000 || options.threads < 1 ||
        options.threads > 256) {
        throw std::invalid_argument(
            "height-gap output requires state/certificate, valid selection and shells, top modes 1..1000, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveHeightGapOutputCli::print_help(std::ostream& out) {
    out << "Projective height-gap shared-output options:\n"
        << "  --state PATH         Fourier state TSV\n"
        << "  --certificate PATH   write English JSON ledger\n"
        << "  --selection NAME     local SLD triad selection\n"
        << "  --first-shell I      lower primitive-height shell\n"
        << "  --second-shell J     upper primitive-height shell\n"
        << "  --top-modes N        retained conjugate output pairs\n"
        << "  --threads N          parallel direct workers\n";
}

int LocalSldProjectiveHeightGapOutputCli::run(
    const LocalSldProjectiveHeightGapOutputCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    const auto report = LocalSldProjectiveHeightGapOutputLedger::analyze(
        state, LocalSldTriadSelection::parse(options.selection),
        options.first_shell, options.second_shell,
        options.top_modes, options.threads);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write height-gap shared-output certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-gap-output-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"first_shell\": " << report.first_shell << ",\n"
        << "  \"second_shell\": " << report.second_shell << ",\n"
        << "  \"shell_gap\": " << report.shell_gap << ",\n"
        << "  \"first_shape_count\": " << report.first_shape_count << ",\n"
        << "  \"second_shape_count\": " << report.second_shape_count << ",\n"
        << "  \"first_interaction_count\": "
        << report.first_interaction_count << ",\n"
        << "  \"second_interaction_count\": "
        << report.second_interaction_count << ",\n"
        << "  \"shared_conjugate_pair_count\": "
        << report.shared_conjugate_pair_count << ",\n"
        << "  \"first_h2_norm2\": "
        << static_cast<double>(report.first_h2_norm2) << ",\n"
        << "  \"second_h2_norm2\": "
        << static_cast<double>(report.second_h2_norm2) << ",\n"
        << "  \"gram_pairing\": "
        << static_cast<double>(report.gram_pairing) << ",\n"
        << "  \"reconstructed_gram_pairing\": "
        << static_cast<double>(report.reconstructed_gram_pairing) << ",\n"
        << "  \"absolute_modal_pairing\": "
        << static_cast<double>(report.absolute_modal_pairing) << ",\n"
        << "  \"correlation\": "
        << static_cast<double>(report.correlation) << ",\n"
        << "  \"first_shared_energy_fraction\": "
        << static_cast<double>(report.first_shared_energy_fraction) << ",\n"
        << "  \"second_shared_energy_fraction\": "
        << static_cast<double>(report.second_shared_energy_fraction) << ",\n"
        << "  \"effective_shared_mode_count\": "
        << static_cast<double>(report.effective_shared_mode_count) << ",\n"
        << "  \"top_signed_gram_fraction\": "
        << static_cast<double>(report.top_signed_gram_fraction) << ",\n"
        << "  \"top_absolute_pairing_fraction\": "
        << static_cast<double>(report.top_absolute_pairing_fraction) << ",\n"
        << "  \"gram_reconstruction_error\": "
        << static_cast<double>(report.gram_reconstruction_error) << ",\n"
        << "  \"exact_gram_reconstruction\": "
        << (report.exact_gram_reconstruction ? "true" : "false") << ",\n"
        << "  \"top_output_modes\": [\n";
    for (std::size_t index = 0; index < report.top_modes.size(); ++index) {
        const auto& mode = report.top_modes[index];
        certificate << "    {\"wave\": [" << mode.wave.x << ", "
            << mode.wave.y << ", " << mode.wave.z << ']'
            << ", \"first_h2_norm2\": "
            << static_cast<double>(mode.first_h2_norm2)
            << ", \"second_h2_norm2\": "
            << static_cast<double>(mode.second_h2_norm2)
            << ", \"pairing_real\": "
            << static_cast<double>(mode.pairing_real)
            << ", \"pairing_imaginary\": "
            << static_cast<double>(mode.pairing_imaginary)
            << ", \"pairing_magnitude\": "
            << static_cast<double>(mode.pairing_magnitude)
            << ", \"modal_correlation\": "
            << static_cast<double>(mode.modal_correlation)
            << ", \"first_energy_fraction\": "
            << static_cast<double>(mode.first_energy_fraction)
            << ", \"second_energy_fraction\": "
            << static_cast<double>(mode.second_energy_fraction)
            << ", \"signed_gram_fraction\": "
            << static_cast<double>(mode.signed_gram_fraction)
            << ", \"absolute_pairing_fraction\": "
            << static_cast<double>(mode.absolute_pairing_fraction)
            << '}'
            << (index + 1 == report.top_modes.size() ? "\n" : ",\n");
    }
    certificate
        << "  ],\n"
        << "  \"finite\": " << (report.finite ? "true" : "false") << ",\n"
        << "  \"finite_ledger_is_not_a_proof\": true,\n"
        << "  \"interpretation\": \"mode-resolved diagnostic of the two shell outputs used by the rejected standalone PNT-13 estimate\"\n"
        << "}\n";
    out << std::setprecision(12)
        << "height-gap shared output shells=" << report.first_shell
        << ',' << report.second_shell
        << " correlation=" << static_cast<double>(report.correlation)
        << " effective_modes="
        << static_cast<double>(report.effective_shared_mode_count)
        << " top_abs_fraction="
        << static_cast<double>(report.top_absolute_pairing_fraction)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.finite ? 0 : 2;
}

}  // namespace lemma
