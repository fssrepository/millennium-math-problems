#include "local_sld_projective_height_gap_triad_attribution.hpp"

#include "local_sld_triad_selection.hpp"
#include "projective_advection_decomposition.hpp"
#include "projective_height_shell_partition.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal vector_norm2(const ComplexVector& value) {
    return std::real(dot_hermitian(value, value));
}

void add(ComplexVector& target, const ComplexVector& source) {
    for (std::size_t component = 0; component < 3; ++component) {
        target[component] += source[component];
    }
}

SpectralReal relative_vector_error(
    const ComplexVector& computed,
    const ComplexVector& expected) {
    ComplexVector difference = computed;
    for (std::size_t component = 0; component < 3; ++component) {
        difference[component] -= expected[component];
    }
    return std::sqrt(std::max(vector_norm2(difference), 0.0L)) /
        std::max(
            std::sqrt(std::max(vector_norm2(expected), 0.0L)),
            1e-30L);
}

ComplexVector interaction_output(
    const SpectralState& state,
    InteractionIndex interaction) {
    const auto [p, q, target] = interaction;
    const SpectralComplex coefficient = SpectralComplex(0.0L, 1.0L) *
        wave_dot(state.waves[q], state.velocity[p]);
    ComplexVector result{};
    for (std::size_t component = 0; component < 3; ++component) {
        result[component] = coefficient * state.velocity[q][component];
    }
    result = project_divergence_free(state.waves[target], result);
    const SpectralReal laplacian = static_cast<SpectralReal>(
        norm_squared(state.waves[target]));
    for (SpectralComplex& component : result) {
        component *= laplacian;
    }
    return result;
}

LocalSldProjectiveHeightGapTriadShellReport analyze_shell(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    const ProjectiveHeightShellGroup& shell,
    std::size_t output_index,
    std::size_t top_interactions,
    int threads) {
    LocalSldProjectiveHeightGapTriadShellReport report;
    report.shell = shell.shell;
    report.minimum_height = shell.minimum_height;
    report.maximum_height = shell.maximum_height;
    report.shape_count = shell.group_indices.size();
    report.shell_interaction_count = shell.interaction_count;
    const SpectralIncrement b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, shell.group_indices,
            state.velocity, state.velocity, threads);
    ComplexVector expected = b[output_index];
    const SpectralReal laplacian = static_cast<SpectralReal>(
        norm_squared(state.waves[output_index]));
    for (SpectralComplex& component : expected) {
        component *= laplacian;
    }
    report.output_norm2 = vector_norm2(expected);

    ComplexVector reconstructed{};
    SpectralReal fraction_square_sum = 0.0L;
    for (const std::size_t group_index : shell.group_indices) {
        const auto& group = groups[group_index];
        for (const InteractionIndex interaction : group.interactions) {
            const auto [p, q, target] = interaction;
            if (target != output_index) {
                continue;
            }
            ++report.target_interaction_count;
            const ComplexVector contribution = interaction_output(
                state, interaction);
            add(reconstructed, contribution);
            LocalSldProjectiveHeightGapTriadContribution row;
            row.primitive_shape = group.primitive_squared_lengths;
            row.advecting_wave = state.waves[p];
            row.advected_wave = state.waves[q];
            row.advecting_energy = vector_norm2(state.velocity[p]);
            row.advected_energy = vector_norm2(state.velocity[q]);
            row.coefficient_magnitude = std::abs(
                wave_dot(state.waves[q], state.velocity[p]));
            row.contribution_norm2 = vector_norm2(contribution);
            row.output_pairing = std::real(dot_hermitian(
                contribution, expected));
            if (report.output_norm2 > 0.0L) {
                row.signed_output_fraction =
                    row.output_pairing / report.output_norm2;
                row.absolute_output_fraction =
                    std::abs(row.output_pairing) / report.output_norm2;
            }
            report.signed_output_fraction_sum +=
                row.signed_output_fraction;
            report.absolute_output_fraction_sum +=
                row.absolute_output_fraction;
            fraction_square_sum += row.absolute_output_fraction *
                row.absolute_output_fraction;
            report.top_contributions.push_back(std::move(row));
        }
    }
    report.reconstructed_output_norm2 = vector_norm2(reconstructed);
    report.output_reconstruction_error = relative_vector_error(
        reconstructed, expected);
    report.exact_output_reconstruction =
        report.output_reconstruction_error < 2e-12L;
    if (fraction_square_sum > 0.0L) {
        report.effective_target_interaction_count =
            report.absolute_output_fraction_sum *
            report.absolute_output_fraction_sum / fraction_square_sum;
    }
    std::sort(
        report.top_contributions.begin(),
        report.top_contributions.end(),
        [](const auto& left, const auto& right) {
            return left.absolute_output_fraction >
                right.absolute_output_fraction;
        });
    if (report.top_contributions.size() > top_interactions) {
        report.top_contributions.resize(top_interactions);
    }
    SpectralReal retained_absolute_fraction = 0.0L;
    for (const auto& row : report.top_contributions) {
        retained_absolute_fraction += row.absolute_output_fraction;
    }
    if (report.absolute_output_fraction_sum > 0.0L) {
        report.top_absolute_fraction_share =
            retained_absolute_fraction /
            report.absolute_output_fraction_sum;
    }
    return report;
}

void write_shell(
    const char* name,
    const LocalSldProjectiveHeightGapTriadShellReport& report,
    std::ostream& output,
    bool trailing_comma) {
    output << std::setprecision(18)
        << "  \"" << name << "\": {\n"
        << "    \"shell\": " << report.shell << ",\n"
        << "    \"height_range\": [" << report.minimum_height
        << ", " << report.maximum_height << "],\n"
        << "    \"shape_count\": " << report.shape_count << ",\n"
        << "    \"shell_interaction_count\": "
        << report.shell_interaction_count << ",\n"
        << "    \"target_interaction_count\": "
        << report.target_interaction_count << ",\n"
        << "    \"output_norm2\": "
        << static_cast<double>(report.output_norm2) << ",\n"
        << "    \"reconstructed_output_norm2\": "
        << static_cast<double>(report.reconstructed_output_norm2) << ",\n"
        << "    \"output_reconstruction_error\": "
        << static_cast<double>(report.output_reconstruction_error) << ",\n"
        << "    \"signed_output_fraction_sum\": "
        << static_cast<double>(report.signed_output_fraction_sum) << ",\n"
        << "    \"absolute_output_fraction_sum\": "
        << static_cast<double>(report.absolute_output_fraction_sum) << ",\n"
        << "    \"effective_target_interaction_count\": "
        << static_cast<double>(report.effective_target_interaction_count)
        << ",\n"
        << "    \"top_absolute_fraction_share\": "
        << static_cast<double>(report.top_absolute_fraction_share) << ",\n"
        << "    \"exact_output_reconstruction\": "
        << (report.exact_output_reconstruction ? "true" : "false")
        << ",\n"
        << "    \"top_interactions\": [\n";
    for (std::size_t index = 0;
         index < report.top_contributions.size(); ++index) {
        const auto& row = report.top_contributions[index];
        output << "      {\"primitive_shape\": ["
            << row.primitive_shape[0] << ", "
            << row.primitive_shape[1] << ", "
            << row.primitive_shape[2] << ']'
            << ", \"advecting_wave\": [" << row.advecting_wave.x
            << ", " << row.advecting_wave.y << ", "
            << row.advecting_wave.z << ']'
            << ", \"advected_wave\": [" << row.advected_wave.x
            << ", " << row.advected_wave.y << ", "
            << row.advected_wave.z << ']'
            << ", \"advecting_energy\": "
            << static_cast<double>(row.advecting_energy)
            << ", \"advected_energy\": "
            << static_cast<double>(row.advected_energy)
            << ", \"coefficient_magnitude\": "
            << static_cast<double>(row.coefficient_magnitude)
            << ", \"contribution_norm2\": "
            << static_cast<double>(row.contribution_norm2)
            << ", \"output_pairing\": "
            << static_cast<double>(row.output_pairing)
            << ", \"signed_output_fraction\": "
            << static_cast<double>(row.signed_output_fraction)
            << ", \"absolute_output_fraction\": "
            << static_cast<double>(row.absolute_output_fraction)
            << '}'
            << (index + 1 == report.top_contributions.size()
                    ? "\n" : ",\n");
    }
    output << "    ]\n"
        << "  }" << (trailing_comma ? ",\n" : "\n");
}

}  // namespace

LocalSldProjectiveHeightGapTriadAttributionReport
LocalSldProjectiveHeightGapTriadAttribution::analyze(
    const SpectralState& state,
    TriadSelection selection,
    int first_shell,
    int second_shell,
    WaveVector output_wave,
    std::size_t top_interactions,
    int threads) {
    if (first_shell < 0 || second_shell <= first_shell ||
        top_interactions < 1 || top_interactions > 1000 ||
        threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "height-gap triad attribution has invalid shells or limits");
    }
    const auto output = state.index.find(output_wave);
    if (output == state.index.end()) {
        throw std::invalid_argument(
            "height-gap triad attribution output wave is absent");
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
            "height-gap triad attribution shell is absent");
    }
    LocalSldProjectiveHeightGapTriadAttributionReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.output_wave = output_wave;
    report.first = analyze_shell(
        state, groups, partition[static_cast<std::size_t>(first_shell)],
        output->second, top_interactions, threads);
    report.second = analyze_shell(
        state, groups, partition[static_cast<std::size_t>(second_shell)],
        output->second, top_interactions, threads);
    report.finite = report.first.exact_output_reconstruction &&
        report.second.exact_output_reconstruction &&
        std::isfinite(report.first.effective_target_interaction_count) &&
        std::isfinite(report.second.effective_target_interaction_count);
    return report;
}

LocalSldProjectiveHeightGapTriadAttributionCliOptions
LocalSldProjectiveHeightGapTriadAttributionCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightGapTriadAttributionCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    bool has_x = false;
    bool has_y = false;
    bool has_z = false;
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
        } else if (name == "--output-x") {
            options.output_wave.x = std::stoi(next(index, name));
            has_x = true;
        } else if (name == "--output-y") {
            options.output_wave.y = std::stoi(next(index, name));
            has_y = true;
        } else if (name == "--output-z") {
            options.output_wave.z = std::stoi(next(index, name));
            has_z = true;
        } else if (name == "--top-interactions") {
            options.top_interactions = static_cast<std::size_t>(
                std::stoull(next(index, name)));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown height-gap triad-attribution option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.first_shell < 0 ||
        options.second_shell <= options.first_shell ||
        options.second_shell > 62 || !has_x || !has_y || !has_z ||
        options.top_interactions < 1 || options.top_interactions > 1000 ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "height-gap triad attribution requires state/certificate, selection, shells, all output coordinates, top interactions 1..1000, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveHeightGapTriadAttributionCli::print_help(
    std::ostream& out) {
    out << "Projective height-gap triad-attribution options:\n"
        << "  --state PATH           Fourier state TSV\n"
        << "  --certificate PATH     write English JSON ledger\n"
        << "  --selection NAME       local SLD triad selection\n"
        << "  --first-shell I        lower primitive-height shell\n"
        << "  --second-shell J       upper primitive-height shell\n"
        << "  --output-x K           target wave x coordinate\n"
        << "  --output-y K           target wave y coordinate\n"
        << "  --output-z K           target wave z coordinate\n"
        << "  --top-interactions N   retained interactions per shell\n"
        << "  --threads N            direct workers\n";
}

int LocalSldProjectiveHeightGapTriadAttributionCli::run(
    const LocalSldProjectiveHeightGapTriadAttributionCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    const auto report = LocalSldProjectiveHeightGapTriadAttribution::analyze(
        state, LocalSldTriadSelection::parse(options.selection),
        options.first_shell, options.second_shell, options.output_wave,
        options.top_interactions, options.threads);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write height-gap triad-attribution certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-gap-triad-attribution-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"output_wave\": [" << report.output_wave.x << ", "
        << report.output_wave.y << ", " << report.output_wave.z << "],\n";
    write_shell("first_shell", report.first, certificate, true);
    write_shell("second_shell", report.second, certificate, true);
    certificate
        << "  \"finite\": " << (report.finite ? "true" : "false") << ",\n"
        << "  \"finite_ledger_is_not_a_proof\": true,\n"
        << "  \"interpretation\": \"ordered-triad attribution of one dominant shared shell output\"\n"
        << "}\n";
    out << std::setprecision(12)
        << "height-gap triad attribution wave=(" << report.output_wave.x
        << ',' << report.output_wave.y << ',' << report.output_wave.z
        << ") effective_interactions="
        << static_cast<double>(
               report.first.effective_target_interaction_count)
        << ',' << static_cast<double>(
               report.second.effective_target_interaction_count)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.finite ? 0 : 2;
}

}  // namespace lemma
