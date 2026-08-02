#include "local_sld_projective_normalization_satellite_scan.hpp"

#include "local_sld_projective_normalization_cauchy_objective.hpp"
#include "local_sld_projective_normalization_objective.hpp"
#include "local_sld_triad_selection.hpp"
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
#include <random>
#include <sstream>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal fitted_slope(
    const std::vector<LocalSldProjectiveNormalizationSatelliteRow>& rows,
    bool actual) {
    SpectralReal n = 0.0L;
    SpectralReal sx = 0.0L;
    SpectralReal sy = 0.0L;
    SpectralReal sxx = 0.0L;
    SpectralReal sxy = 0.0L;
    for (const auto& row : rows) {
        const SpectralReal value = actual
            ? row.selected_channel_power_one
            : row.cauchy_bound_power_one;
        if (row.cutoff < 1 || !(value > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
        const SpectralReal y = std::log(value);
        n += 1.0L;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    const SpectralReal denominator = n * sxx - sx * sx;
    return n >= 2.0L && std::abs(denominator) > 1e-30L
        ? (n * sxy - sx * sy) / denominator
        : 0.0L;
}

SpectralState outer_shell_satellite(
    int cutoff,
    std::uint64_t seed,
    int mode_pairs) {
    std::mt19937_64 generator(seed);
    SpectralState state = SpectralStateFactory::random(cutoff, generator);
    if (mode_pairs == 0) {
        const int minimum_shell = std::max(1, (cutoff + 2) / 3);
        const int maximum_shell = std::max(
            minimum_shell, 2 * cutoff / 3);
        for (std::size_t index = 0; index < state.waves.size(); ++index) {
            const WaveVector wave = state.waves[index];
            const int radius = std::max({
                std::abs(wave.x), std::abs(wave.y), std::abs(wave.z)});
            if (radius < minimum_shell || radius > maximum_shell) {
                state.velocity[index] = {};
            }
        }
        SpectralStateOps::normalize_energy(state);
        return state;
    }
    for (auto& velocity : state.velocity) {
        velocity = {};
    }
    const int high = cutoff - 1;
    const std::vector<WaveVector> structured = {
        {high, 0, 0}, {high, 1, 0}, {high, 0, 1}, {high, 1, 1},
        {0, high, 0}, {1, high, 0}, {0, high, 1}, {1, high, 1},
        {0, 0, high}, {1, 0, high}, {0, 1, high}, {1, 1, high},
        {high, -1, 0}, {high, 0, -1}, {0, high, -1}, {0, -1, high},
    };
    if (structured.size() < static_cast<std::size_t>(mode_pairs)) {
        throw std::runtime_error("not enough outer-shell satellite modes");
    }
    for (int pair = 0; pair < mode_pairs; ++pair) {
        const WaveVector wave = structured[static_cast<std::size_t>(pair)];
        ComplexVector value = {
            SpectralComplex{1.0L, 0.5L},
            SpectralComplex{-0.25L, 1.0L},
            SpectralComplex{0.75L, -0.5L},
        };
        value = project_divergence_free(wave, value);
        const std::size_t index = state.index.at(wave);
        state.velocity[index] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

void write_json(
    const LocalSldProjectiveNormalizationSatelliteReport& report,
    const LocalSldProjectiveNormalizationSatelliteOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-low-high-satellite-v1\",\n"
        << "  \"base_state_path\": \"" << options.base_state_path
        << "\",\n"
        << "  \"base_cutoff\": " << options.base_cutoff << ",\n"
        << "  \"core_maximum_height\": "
        << options.core_maximum_height << ",\n"
        << "  \"satellite_amplitude_formula\": \"coefficient/cutoff^2\",\n"
        << "  \"satellite_layout\": \""
        << (options.satellite_mode_pairs == 0
                ? "dense random middle-frequency band"
                : "structured scaled outer-shell wave family")
        << "\",\n"
        << "  \"satellite_coefficient\": "
        << static_cast<double>(options.satellite_coefficient) << ",\n"
        << "  \"satellite_mode_pairs\": "
        << options.satellite_mode_pairs << ",\n"
        << "  \"seed\": " << options.seed << ",\n"
        << "  \"triad_selection\": \"" << options.selection << "\",\n"
        << "  \"fitted_cauchy_cutoff_slope\": "
        << static_cast<double>(report.fitted_cauchy_cutoff_slope) << ",\n"
        << "  \"fitted_actual_cutoff_slope\": "
        << static_cast<double>(report.fitted_actual_cutoff_slope) << ",\n"
        << "  \"maximum_cauchy_bound\": "
        << static_cast<double>(report.maximum_cauchy_bound) << ",\n"
        << "  \"maximum_selected_channel\": "
        << static_cast<double>(report.maximum_selected_channel) << ",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        output << "    {\"cutoff\": " << row.cutoff
            << ", \"satellite_amplitude\": "
            << static_cast<double>(row.satellite_amplitude)
            << ", \"satellite_energy_fraction\": "
            << static_cast<double>(row.satellite_energy_fraction)
            << ", \"enstrophy\": "
            << static_cast<double>(row.enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(row.palinstrophy)
            << ", \"full_stretching\": "
            << static_cast<double>(row.full_stretching)
            << ", \"selected_aggregate_h1_norm2\": "
            << static_cast<double>(row.selected_aggregate_h1_norm2)
            << ", \"tail_aggregate_h2_norm2\": "
            << static_cast<double>(row.tail_aggregate_h2_norm2)
            << ", \"cauchy_bound_power_one\": "
            << static_cast<double>(row.cauchy_bound_power_one)
            << ", \"selected_channel_power_one\": "
            << static_cast<double>(row.selected_channel_power_one)
            << ", \"actual_over_cauchy_bound\": "
            << static_cast<double>(row.actual_over_cauchy_bound)
            << ", \"quarter_compensated_bound\": "
            << static_cast<double>(row.quarter_compensated_bound)
            << ", \"state_path\": \"" << row.state_path
            << "\", \"finite\": " << (row.finite ? "true" : "false")
            << '}' << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_row_finite\": "
        << (report.every_row_finite ? "true" : "false") << ",\n"
        << "  \"finite_satellite_scan_is_not_a_proof\": true,\n"
        << "  \"uniform_cauchy_majorant_proved\": false,\n"
        << "  \"interpretation\": \"finite locality stress test only: structured scale-separated satellites may be excluded exactly by the local triad selector, while a dense comparable-frequency band tests local high-high output; neither outcome proves or falsifies a cutoff-uniform majorant without an unbounded-cutoff construction\"\n"
        << "}\n";
}

}  // namespace

LocalSldProjectiveNormalizationSatelliteOptions
LocalSldProjectiveNormalizationSatelliteCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationSatelliteOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--base-state") {
            options.base_state_path = next(index, name);
        } else if (name == "--base-cutoff") {
            options.base_cutoff = std::stoi(next(index, name));
        } else if (name == "--cutoff") {
            options.cutoffs.push_back(std::stoi(next(index, name)));
        } else if (name == "--projective-core-height") {
            options.core_maximum_height = static_cast<SpectralInteger>(
                std::stoll(next(index, name)));
        } else if (name == "--satellite-coefficient") {
            options.satellite_coefficient = std::stold(next(index, name));
        } else if (name == "--satellite-modes") {
            options.satellite_mode_pairs = std::stoi(next(index, name));
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--state-dir") {
            options.state_directory = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown normalization-satellite option: " + name);
        }
    }
    std::sort(options.cutoffs.begin(), options.cutoffs.end());
    options.cutoffs.erase(
        std::unique(options.cutoffs.begin(), options.cutoffs.end()),
        options.cutoffs.end());
    if (options.base_state_path.empty() ||
        options.state_directory.empty() ||
        options.certificate_path.empty() ||
        options.cutoffs.size() < 2 ||
        options.base_cutoff < 1 || options.base_cutoff > 10 ||
        options.core_maximum_height < 1 ||
        !(options.satellite_coefficient > 0.0L) ||
        !std::isfinite(options.satellite_coefficient) ||
        options.satellite_mode_pairs < 0 ||
        options.satellite_mode_pairs > 64 ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.threads < 1 || options.threads > 256 ||
        std::any_of(
            options.cutoffs.begin(), options.cutoffs.end(),
            [&](int cutoff) {
                return cutoff <= options.base_cutoff + 1 || cutoff > 12 ||
                    options.satellite_coefficient /
                        static_cast<SpectralReal>(cutoff * cutoff) >= 1.0L;
            })) {
        throw std::invalid_argument(
            "normalization-satellite requires a base state, output paths, at least two cutoffs above base+1 through 12, positive core height/coefficient, valid selection, and threads 1..256");
    }
    return options;
}

void LocalSldProjectiveNormalizationSatelliteCli::print_help(
    std::ostream& out) {
    out << "Projective normalization low/high satellite scan options:\n"
        << "  --base-state PATH           source state projected to the low core\n"
        << "  --base-cutoff K             low-frequency projection cutoff\n"
        << "  --cutoff K                  target cutoff; repeatable\n"
        << "  --projective-core-height H  fixed primitive-height core\n"
        << "  --satellite-coefficient C   outer-shell amplitude C/K^2\n"
        << "  --satellite-modes N         1..64 structured pairs; 0 selects a dense middle band\n"
        << "  --selection NAME            local SLD triad selection\n"
        << "  --seed N                    deterministic satellite seed\n"
        << "  --threads N                 direct workers\n"
        << "  --state-dir PATH            write generated states\n"
        << "  --certificate PATH          write English JSON scan\n";
}

int LocalSldProjectiveNormalizationSatelliteCli::run(
    const LocalSldProjectiveNormalizationSatelliteOptions& options,
    std::ostream& out) {
    SpectralState base = SpectralStateFactory::project(
        SpectralStateReader::read_tsv(options.base_state_path),
        options.base_cutoff);
    SpectralStateOps::normalize_energy(base);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const TriadSelection selection = LocalSldTriadSelection::parse(
        options.selection);
    LocalSldProjectiveNormalizationSatelliteReport report;
    report.every_row_finite = true;
    std::filesystem::create_directories(options.state_directory);
    for (const int cutoff : options.cutoffs) {
        std::mt19937_64 layout_generator(options.seed);
        const SpectralState lifted_base = SpectralStateFactory::lift(
            base, cutoff, layout_generator);
        const SpectralState satellite = outer_shell_satellite(
            cutoff, options.seed + static_cast<std::uint64_t>(cutoff),
            options.satellite_mode_pairs);
        const SpectralReal amplitude = options.satellite_coefficient /
            static_cast<SpectralReal>(cutoff * cutoff);
        const SpectralState state = SpectralStateBlend::blend_on_energy_sphere(
            lifted_base, satellite, amplitude);
        const auto cauchy = LocalSldProjectiveNormalizationCauchyObjective(
            dynamics, selection, options.core_maximum_height,
            options.threads).evaluate(state);
        const auto actual = LocalSldProjectiveNormalizationObjective(
            dynamics, selection, options.core_maximum_height,
            options.threads,
            LocalSldProjectiveNormalizationComponent::
                selected_stretching_tail_cross).evaluate(state);
        LocalSldProjectiveNormalizationSatelliteRow row;
        row.cutoff = cutoff;
        row.satellite_amplitude = amplitude;
        row.satellite_energy_fraction = amplitude * amplitude;
        row.enstrophy = cauchy.enstrophy;
        row.palinstrophy = cauchy.palinstrophy;
        row.full_stretching = cauchy.full_stretching;
        row.selected_aggregate_h1_norm2 =
            cauchy.selected_aggregate_h1_norm2;
        row.tail_aggregate_h2_norm2 = cauchy.tail_aggregate_h2_norm2;
        row.cauchy_bound_power_one = cauchy.cauchy_bound_power_one;
        row.selected_channel_power_one =
            actual.selected_stretching_tail_cross_power_one;
        if (row.cauchy_bound_power_one > 0.0L) {
            row.actual_over_cauchy_bound =
                row.selected_channel_power_one /
                row.cauchy_bound_power_one;
        }
        row.quarter_compensated_bound = std::pow(
            static_cast<SpectralReal>(options.core_maximum_height), 0.25L) *
            row.cauchy_bound_power_one;
        row.state_path = (std::filesystem::path(options.state_directory) /
            ("K" + std::to_string(cutoff) + ".tsv")).string();
        row.finite = cauchy.finite && actual.finite &&
            std::isfinite(row.actual_over_cauchy_bound);
        std::ostringstream metadata;
        metadata << std::setprecision(18)
            << "projective normalization low/high satellite; base="
            << options.base_state_path
            << "; base_cutoff=" << options.base_cutoff
            << "; cutoff=" << cutoff
            << "; satellite_amplitude="
            << static_cast<double>(amplitude)
            << "; satellite_mode_pairs="
            << options.satellite_mode_pairs
            << "; candidate_lemma_proved=false";
        SpectralStateWriter::write_tsv(
            row.state_path, state, metadata.str());
        report.every_row_finite = report.every_row_finite && row.finite;
        report.maximum_cauchy_bound = std::max(
            report.maximum_cauchy_bound, row.cauchy_bound_power_one);
        report.maximum_selected_channel = std::max(
            report.maximum_selected_channel,
            row.selected_channel_power_one);
        report.rows.push_back(std::move(row));
    }
    report.fitted_cauchy_cutoff_slope = fitted_slope(report.rows, false);
    report.fitted_actual_cutoff_slope = fitted_slope(report.rows, true);
    const std::filesystem::path certificate_path(options.certificate_path);
    if (!certificate_path.parent_path().empty()) {
        std::filesystem::create_directories(certificate_path.parent_path());
    }
    std::ofstream certificate(certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write normalization-satellite certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective normalization satellite scan rows="
        << report.rows.size()
        << " cauchy_slope="
        << static_cast<double>(report.fitted_cauchy_cutoff_slope)
        << " actual_slope="
        << static_cast<double>(report.fitted_actual_cutoff_slope)
        << " max_cauchy="
        << static_cast<double>(report.maximum_cauchy_bound)
        << " max_actual="
        << static_cast<double>(report.maximum_selected_channel) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.every_row_finite ? 0 : 2;
}

}  // namespace lemma
