#include "local_sld_projective_height_coercivity_path_scan.hpp"

#include "local_quartic_closure_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "projective_height_envelope_kernel.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_state_blend.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

SpectralReal tail_slope(
    const std::vector<LocalSldProjectiveHeightCoercivityPathRow>& rows,
    SpectralReal LocalSldProjectiveHeightCoercivityPathRow::*field) {
    const std::size_t count = std::max<std::size_t>(
        3, (rows.size() + 1) / 2);
    const std::size_t first = rows.size() - std::min(count, rows.size());
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    std::size_t valid = 0;
    for (std::size_t index = first; index < rows.size(); ++index) {
        const SpectralReal x_value = rows[index].epsilon;
        const SpectralReal y_value = rows[index].*field;
        if (!(x_value > 0.0L) || !(y_value > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(x_value);
        const SpectralReal y = std::log(y_value);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        ++valid;
    }
    if (valid < 2) {
        return std::numeric_limits<SpectralReal>::quiet_NaN();
    }
    const SpectralReal n = static_cast<SpectralReal>(valid);
    const SpectralReal denominator = n * sum_xx - sum_x * sum_x;
    return std::abs(denominator) > 1e-30L
        ? (n * sum_xy - sum_x * sum_y) / denominator
        : std::numeric_limits<SpectralReal>::quiet_NaN();
}

WaveVector parse_wave(std::string value) {
    std::replace(value.begin(), value.end(), ',', ' ');
    std::istringstream input(value);
    WaveVector wave;
    if (!(input >> wave.x >> wave.y >> wave.z) ||
        norm_squared(wave) == 0) {
        throw std::invalid_argument(
            "base wave must be a nonzero comma-separated integer triple");
    }
    std::string trailing;
    if (input >> trailing) {
        throw std::invalid_argument("base wave has trailing data");
    }
    return wave;
}

void write_json(
    const LocalSldProjectiveHeightCoercivityPathReport& report,
    const LocalSldProjectiveHeightCoercivityPathScanOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-coercivity-path-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"base_wave\": [" << report.base_wave.x << ", "
        << report.base_wave.y << ", " << report.base_wave.z << "],\n"
        << "  \"base_selected_advection_h1_weight\": "
        << static_cast<double>(report.base_selected_advection_h1_weight)
        << ",\n"
        << "  \"base_is_selected_advection_null\": "
        << (report.base_is_selected_advection_null ? "true" : "false")
        << ",\n"
        << "  \"height_shell_count\": " << report.shell_count << ",\n"
        << "  \"tail_log_log_slopes\": {\n"
        << "    \"paired_envelope\": "
        << static_cast<double>(report.paired_envelope_tail_slope) << ",\n"
        << "    \"outer_h1_weight\": "
        << static_cast<double>(report.outer_h1_weight_tail_slope) << ",\n"
        << "    \"coercivity_ratio\": "
        << static_cast<double>(report.ratio_tail_slope) << "\n"
        << "  },\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"epsilon\": " << static_cast<double>(row.epsilon)
            << ", \"energy\": " << static_cast<double>(row.energy)
            << ", \"paired_envelope\": "
            << static_cast<double>(row.paired_envelope)
            << ", \"outer_h1_weight\": "
            << static_cast<double>(row.outer_h1_weight)
            << ", \"coercivity_ratio\": "
            << static_cast<double>(row.coercivity_ratio)
            << ", \"squared_coercivity_ratio\": "
            << static_cast<double>(row.squared_coercivity_ratio)
            << ", \"signed_bracket\": "
            << static_cast<double>(row.signed_bracket)
            << ", \"component_envelopes\": [";
        for (std::size_t component = 0;
             component < row.component_envelopes.size(); ++component) {
            output << static_cast<double>(
                row.component_envelopes[component]);
            if (component + 1 < row.component_envelopes.size()) {
                output << ", ";
            }
        }
        output << "]}" << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"finite_path_scan_is_not_a_proof\": true,\n"
        << "  \"candidate_outer_coercivity_lemma_proved\": false,\n"
        << "  \"interpretation\": \"a negative ratio slope indicates growth toward the selected-advection null state\"\n"
        << "}\n";
}

}  // namespace

SpectralState
LocalSldProjectiveHeightCoercivityPathScan::monochromatic_base(
    const SpectralState& layout,
    WaveVector wave) {
    if (norm_squared(wave) == 0) {
        throw std::invalid_argument("monochromatic base wave is zero");
    }
    const auto positive = layout.index.find(wave);
    const auto negative = layout.index.find(-wave);
    if (positive == layout.index.end() || negative == layout.index.end()) {
        throw std::invalid_argument(
            "monochromatic base wave is absent from the state layout");
    }
    SpectralState result = layout;
    for (ComplexVector& velocity : result.velocity) {
        velocity = {};
    }
    ComplexVector value{
        SpectralComplex{1.0L, 0.0L},
        SpectralComplex{0.0L, 1.0L},
        SpectralComplex{0.5L, -0.25L}};
    value = project_divergence_free(wave, value);
    result.velocity[positive->second] = value;
    result.velocity[negative->second] = conjugate(value);
    SpectralStateOps::normalize_energy(result);
    return result;
}

LocalSldProjectiveHeightCoercivityPathReport
LocalSldProjectiveHeightCoercivityPathScan::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& perturbation,
    TriadSelection selection,
    WaveVector base_wave,
    SpectralReal minimum_epsilon,
    SpectralReal maximum_epsilon,
    int samples,
    int threads) {
    if (!(minimum_epsilon > 0.0L) ||
        !(maximum_epsilon >= minimum_epsilon) ||
        !(maximum_epsilon <= 1.0L) || samples < 2 || samples > 1000 ||
        threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "invalid projective-height coercivity path parameters");
    }
    const SpectralState base = monochromatic_base(
        perturbation, base_wave);
    LocalSldProjectiveHeightCoercivityPathReport report;
    report.cutoff = SpectralStateOps::cutoff(perturbation);
    report.base_wave = base_wave;

    const LocalQuarticClosureObjective objective(dynamics, selection);
    const LocalQuarticClosureObjectiveValue base_selected =
        objective.evaluate(base);
    const ProjectiveHeightEnvelopeMoment base_envelope =
        ProjectiveHeightEnvelopeKernel::evaluate(
            base, selection, base_selected.enstrophy,
            base_selected.palinstrophy, false, threads, true);
    report.base_selected_advection_h1_weight = base_envelope.outer_h1_sum;
    report.base_is_selected_advection_null =
        std::abs(base_envelope.outer_h1_sum) < 1e-28L &&
        std::abs(base_envelope.absolute_component_envelope) < 1e-28L;

    report.rows.reserve(static_cast<std::size_t>(samples));
    const SpectralReal log_minimum = std::log(minimum_epsilon);
    const SpectralReal log_maximum = std::log(maximum_epsilon);
    for (int index = 0; index < samples; ++index) {
        const SpectralReal fraction = static_cast<SpectralReal>(index) /
            static_cast<SpectralReal>(samples - 1);
        const SpectralReal epsilon = std::exp(
            log_maximum + fraction * (log_minimum - log_maximum));
        const SpectralState state = SpectralStateBlend::blend_on_energy_sphere(
            base, perturbation, epsilon);
        const LocalQuarticClosureObjectiveValue selected =
            objective.evaluate(state);
        const ProjectiveHeightEnvelopeMoment envelope =
            ProjectiveHeightEnvelopeKernel::evaluate(
                state, selection, selected.enstrophy,
                selected.palinstrophy, false, threads, true);
        LocalSldProjectiveHeightCoercivityPathRow row;
        row.epsilon = epsilon;
        row.energy = SpectralStateOps::energy(state);
        row.paired_envelope = envelope.absolute_component_envelope;
        row.outer_h1_weight = envelope.outer_h1_sum;
        row.signed_bracket = selected.signed_two_entry_bracket;
        row.component_envelopes = envelope.absolute_component_sums;
        if (row.outer_h1_weight > 1e-30L) {
            row.coercivity_ratio =
                row.paired_envelope / row.outer_h1_weight;
            row.squared_coercivity_ratio =
                row.coercivity_ratio * row.coercivity_ratio;
        }
        report.shell_count = envelope.height_shell_count;
        report.rows.push_back(row);
    }
    report.paired_envelope_tail_slope = tail_slope(
        report.rows,
        &LocalSldProjectiveHeightCoercivityPathRow::paired_envelope);
    report.outer_h1_weight_tail_slope = tail_slope(
        report.rows,
        &LocalSldProjectiveHeightCoercivityPathRow::outer_h1_weight);
    report.ratio_tail_slope = tail_slope(
        report.rows,
        &LocalSldProjectiveHeightCoercivityPathRow::coercivity_ratio);
    return report;
}

LocalSldProjectiveHeightCoercivityPathScanOptions
LocalSldProjectiveHeightCoercivityPathScanCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightCoercivityPathScanOptions options;
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
        } else if (name == "--base-wave") {
            options.base_wave = parse_wave(next(index, name));
        } else if (name == "--min-epsilon") {
            options.minimum_epsilon = std::stold(next(index, name));
        } else if (name == "--max-epsilon") {
            options.maximum_epsilon = std::stold(next(index, name));
        } else if (name == "--samples") {
            options.samples = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown projective-height coercivity-path option: " +
                name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        !(options.minimum_epsilon > 0.0L) ||
        !(options.maximum_epsilon >= options.minimum_epsilon) ||
        !(options.maximum_epsilon <= 1.0L) ||
        options.samples < 2 || options.samples > 1000 ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "coercivity-path requires state, certificate, valid selection, 0 < min <= max <= 1, samples 2..1000, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveHeightCoercivityPathScanCli::print_help(
    std::ostream& out) {
    out << "Projective-height outer-coercivity null-path options:\n"
        << "  --state PATH          perturbing Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON scan\n"
        << "  --selection NAME      local SLD triad selection\n"
        << "  --base-wave X,Y,Z     monochromatic null-state wave\n"
        << "  --min-epsilon X       smallest perturbation amplitude\n"
        << "  --max-epsilon X       largest perturbation amplitude\n"
        << "  --samples N           logarithmic path samples\n"
        << "  --threads N           parallel projective workers\n";
}

int LocalSldProjectiveHeightCoercivityPathScanCli::run(
    const LocalSldProjectiveHeightCoercivityPathScanOptions& options,
    std::ostream& out) {
    const SpectralState perturbation = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const auto report = LocalSldProjectiveHeightCoercivityPathScan::analyze(
        dynamics, perturbation,
        LocalSldTriadSelection::parse(options.selection),
        options.base_wave, options.minimum_epsilon,
        options.maximum_epsilon, options.samples, options.threads);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective-height coercivity-path certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective-height coercivity path cutoff=" << report.cutoff
        << " base_null="
        << (report.base_is_selected_advection_null ? "true" : "false")
        << " paired_slope="
        << static_cast<double>(report.paired_envelope_tail_slope)
        << " outer_slope="
        << static_cast<double>(report.outer_h1_weight_tail_slope)
        << " ratio_slope="
        << static_cast<double>(report.ratio_tail_slope)
        << " terminal_ratio="
        << static_cast<double>(report.rows.back().coercivity_ratio)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.base_is_selected_advection_null ? 0 : 2;
}

}  // namespace lemma
