#include "local_sld_projective_normalization_alternating_adversary.hpp"

#include "gradient_adversary.hpp"
#include "local_sld_projective_normalization_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "spectral_adjoint.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_objective.hpp"
#include "state_analysis.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

std::string local_component_objective(const std::string& component) {
    if (component == "core-stretching-tail-cross") {
        return "local-projective-core-stretching-tail-cross-ratio";
    }
    if (component == "tail-stretching-core-cross") {
        return "local-projective-tail-stretching-core-cross-ratio";
    }
    if (component == "tail-stretching-tail-cross") {
        return "local-projective-tail-stretching-tail-cross-ratio";
    }
    throw std::invalid_argument(
        "unknown normalization alternating component: " + component);
}

void write_json(
    const LocalSldProjectiveNormalizationAlternatingReport& report,
    const LocalSldProjectiveNormalizationAlternatingOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-alternating-adversary-v1\",\n"
        << "  \"initial_state_path\": \"" << options.state_path
        << "\",\n"
        << "  \"output_state_path\": \""
        << options.output_state_path << "\",\n"
        << "  \"triad_selection\": \"" << options.selection
        << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"core_maximum_height\": "
        << report.core_maximum_height << ",\n"
        << "  \"component\": \"" << report.component << "\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"cycles\": " << options.cycles << ",\n"
        << "  \"component_iterations\": "
        << options.component_iterations << ",\n"
        << "  \"open_iterations\": "
        << options.open_iterations << ",\n"
        << "  \"initial_open_power_one\": "
        << static_cast<double>(report.initial_open_power_one) << ",\n"
        << "  \"final_open_power_one\": "
        << static_cast<double>(report.final_open_power_one) << ",\n"
        << "  \"improvement_factor\": "
        << static_cast<double>(report.improvement_factor) << ",\n"
        << "  \"phases\": [\n";
    for (std::size_t index = 0; index < report.phases.size(); ++index) {
        const auto& phase = report.phases[index];
        output << "    {\"cycle\": " << phase.cycle
            << ", \"component_initial\": "
            << static_cast<double>(phase.component_initial)
            << ", \"component_final\": "
            << static_cast<double>(phase.component_final)
            << ", \"open_initial\": "
            << static_cast<double>(phase.open_initial)
            << ", \"open_final\": "
            << static_cast<double>(phase.open_final)
            << ", \"component_projected_gradient_norm\": "
            << static_cast<double>(
                   phase.component_projected_gradient_norm)
            << ", \"open_projected_gradient_norm\": "
            << static_cast<double>(phase.open_projected_gradient_norm)
            << ", \"component_accepted_steps\": "
            << phase.component_accepted_steps
            << ", \"open_accepted_steps\": "
            << phase.open_accepted_steps
            << ", \"component_evaluations\": "
            << phase.component_evaluations
            << ", \"open_evaluations\": "
            << phase.open_evaluations << '}'
            << (index + 1 == report.phases.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"every_phase_finite\": "
        << (report.every_phase_finite ? "true" : "false") << ",\n"
        << "  \"finite_alternating_search_is_not_a_proof\": true,\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"interpretation\": \"alternating an exact component gradient with the exact open-sum gradient exposes maxima missed by a single monolithic L-BFGS trajectory\"\n"
        << "}\n";
}

GradientSearchOptions search_options(
    const LocalSldProjectiveNormalizationAlternatingOptions& options,
    TriadSelection selection,
    const std::string& objective,
    int iterations) {
    GradientSearchOptions search;
    search.iterations = iterations;
    search.line_search_steps = options.line_search_steps;
    search.initial_step = options.initial_step;
    search.method = "lbfgs";
    search.lbfgs_history = options.lbfgs_history;
    search.objective = objective;
    search.closure_selection = selection;
    search.projective_core_maximum_height =
        options.core_maximum_height;
    search.objective_threads = options.threads;
    return search;
}

}  // namespace

LocalSldProjectiveNormalizationAlternatingOptions
LocalSldProjectiveNormalizationAlternatingAdversaryCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationAlternatingOptions options;
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
        } else if (name == "--output-state") {
            options.output_state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--component") {
            options.component = next(index, name);
        } else if (name == "--projective-core-height") {
            options.core_maximum_height =
                static_cast<SpectralInteger>(
                    std::stoll(next(index, name)));
        } else if (name == "--cycles") {
            options.cycles = std::stoi(next(index, name));
        } else if (name == "--component-iterations") {
            options.component_iterations = std::stoi(next(index, name));
        } else if (name == "--open-iterations") {
            options.open_iterations = std::stoi(next(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next(index, name));
        } else if (name == "--lbfgs-history") {
            options.lbfgs_history = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown normalization-alternating option: " + name);
        }
    }
    const bool valid_component =
        options.component == "core-stretching-tail-cross" ||
        options.component == "tail-stretching-core-cross" ||
        options.component == "tail-stretching-tail-cross";
    if (options.state_path.empty() ||
        options.output_state_path.empty() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        !valid_component || options.core_maximum_height < 1 ||
        options.cycles < 1 || options.cycles > 100 ||
        options.component_iterations < 0 ||
        options.open_iterations < 0 ||
        options.line_search_steps < 1 ||
        options.lbfgs_history < 1 || options.lbfgs_history > 64 ||
        options.threads < 1 || options.threads > 256 ||
        !(options.initial_step > 0.0L) ||
        !std::isfinite(options.initial_step)) {
        throw std::invalid_argument(
            "normalization-alternating requires state/output/certificate, valid selection/component, positive core height, cycles 1..100, nonnegative phase iterations, and valid search parameters");
    }
    return options;
}

void LocalSldProjectiveNormalizationAlternatingAdversaryCli::print_help(
    std::ostream& out) {
    out << "Projective normalization alternating exact-gradient options:\n"
        << "  --state PATH                 initial Fourier TSV\n"
        << "  --output-state PATH          final Fourier TSV\n"
        << "  --certificate PATH           write English JSON trace\n"
        << "  --selection NAME             local SLD triad selection\n"
        << "  --component NAME             core-stretching-tail-cross, tail-stretching-core-cross, or tail-stretching-tail-cross\n"
        << "  --projective-core-height H   fixed primitive-height core\n"
        << "  --cycles N                   alternating phase pairs\n"
        << "  --component-iterations N     L-BFGS component steps per cycle\n"
        << "  --open-iterations N          L-BFGS open-sum steps per cycle\n"
        << "  --line-search N              trials per L-BFGS step\n"
        << "  --lbfgs-history N            curvature pairs\n"
        << "  --threads N                  direct workers\n"
        << "  --step X                     initial line-search step\n";
}

int LocalSldProjectiveNormalizationAlternatingAdversaryCli::run(
    const LocalSldProjectiveNormalizationAlternatingOptions& options,
    std::ostream& out) {
    SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const SpectralObjective spectral_objective(dynamics);
    const SpectralAdjoint adjoint(dynamics, spectral_objective);
    const GradientAdversary adversary(
        dynamics, spectral_objective, adjoint);
    const TriadSelection selection =
        LocalSldTriadSelection::parse(options.selection);
    const LocalSldProjectiveNormalizationObjective open_objective(
        dynamics, selection, options.core_maximum_height,
        options.threads);

    LocalSldProjectiveNormalizationAlternatingReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.core_maximum_height = options.core_maximum_height;
    report.component = options.component;
    report.initial_open_power_one = open_objective.evaluate(state)
        .palinstrophy_normalization_power_one;
    report.every_phase_finite = true;
    const std::string component_objective =
        local_component_objective(options.component);
    for (int cycle = 0; cycle < options.cycles; ++cycle) {
        const GradientSearchResult component = adversary.maximize_q(
            state,
            search_options(
                options, selection, component_objective,
                options.component_iterations));
        state = component.state;
        const GradientSearchResult open = adversary.maximize_q(
            state,
            search_options(
                options, selection,
                "local-projective-open-palinstrophy-normalization-ratio",
                options.open_iterations));
        state = open.state;
        LocalSldProjectiveNormalizationAlternatingPhase phase;
        phase.cycle = cycle + 1;
        phase.component_initial = std::sqrt(std::max(
            0.0L, component.initial_objective));
        phase.component_final = std::sqrt(std::max(
            0.0L, component.objective));
        phase.open_initial = std::sqrt(std::max(
            0.0L, open.initial_objective));
        phase.open_final = std::sqrt(std::max(
            0.0L, open.objective));
        phase.component_projected_gradient_norm =
            component.final_projected_gradient_norm;
        phase.open_projected_gradient_norm =
            open.final_projected_gradient_norm;
        phase.component_accepted_steps = component.accepted_steps;
        phase.open_accepted_steps = open.accepted_steps;
        phase.component_evaluations = component.trajectory_evaluations;
        phase.open_evaluations = open.trajectory_evaluations;
        report.every_phase_finite = report.every_phase_finite &&
            std::isfinite(phase.component_final) &&
            std::isfinite(phase.open_final);
        report.phases.push_back(phase);
    }
    report.state = state;
    report.final_open_power_one = open_objective.evaluate(state)
        .palinstrophy_normalization_power_one;
    if (report.initial_open_power_one > 0.0L) {
        report.improvement_factor = report.final_open_power_one /
            report.initial_open_power_one;
    }

    std::ostringstream metadata;
    metadata << std::setprecision(18)
        << "projective normalization alternating exact-gradient winner; component="
        << options.component
        << "; core_maximum_height=" << options.core_maximum_height
        << "; cycles=" << options.cycles
        << "; initial_open="
        << static_cast<double>(report.initial_open_power_one)
        << "; final_open="
        << static_cast<double>(report.final_open_power_one)
        << "; candidate_lemma_proved=false";
    SpectralStateWriter::write_tsv(
        options.output_state_path, state, metadata.str());
    const std::filesystem::path certificate_path(
        options.certificate_path);
    if (!certificate_path.parent_path().empty()) {
        std::filesystem::create_directories(
            certificate_path.parent_path());
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write normalization alternating certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "projective normalization alternating cutoff=" << report.cutoff
        << " component=" << report.component
        << " open="
        << static_cast<double>(report.initial_open_power_one)
        << " -> "
        << static_cast<double>(report.final_open_power_one)
        << " factor="
        << static_cast<double>(report.improvement_factor) << '\n'
        << "Certificate written to " << options.certificate_path << '\n'
        << "State written to " << options.output_state_path << '\n';
    return report.every_phase_finite ? 0 : 2;
}

}  // namespace lemma
