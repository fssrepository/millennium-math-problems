#include "local_sld_projective_height_gap_correlation_adversary.hpp"

#include "gradient_adversary.hpp"
#include "local_sld_projective_height_gap_correlation_objective.hpp"
#include "local_sld_triad_selection.hpp"
#include "spectral_adjoint.hpp"
#include "spectral_galerkin.hpp"
#include "spectral_objective.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

struct GradientDirectionalCheck {
    SpectralReal absolute_error = 0.0L;
    SpectralReal relative_error = 0.0L;
};

GradientDirectionalCheck gradient_directional_check(
    const SpectralDynamics& dynamics,
    const LocalSldProjectiveHeightGapCorrelationObjective& objective,
    const SpectralState& state,
    std::uint64_t seed) {
    std::mt19937_64 generator(seed);
    const SpectralState random = SpectralStateFactory::random(
        SpectralStateOps::cutoff(state), generator);
    SpectralIncrement direction = random.velocity;
    const SpectralReal norm = std::sqrt(std::max(
        pairing(direction, direction), 0.0L));
    if (!(norm > 0.0L)) {
        return {};
    }
    for (ComplexVector& mode : direction) {
        for (SpectralComplex& component : mode) {
            component /= norm;
        }
    }
    const SpectralIncrement gradient = objective.gradient(state);
    const SpectralReal analytic = pairing(gradient, direction);
    const SpectralReal epsilon = 2e-6L;
    SpectralState plus = dynamics.add_increment(
        state, direction, epsilon);
    SpectralState minus = dynamics.add_increment(
        state, direction, -epsilon);
    dynamics.enforce_constraints(plus);
    dynamics.enforce_constraints(minus);
    const SpectralReal central =
        (objective.evaluate(plus).weighted_correlation_squared -
         objective.evaluate(minus).weighted_correlation_squared) /
        (2.0L * epsilon);
    GradientDirectionalCheck result;
    result.absolute_error = std::abs(central - analytic);
    result.relative_error = result.absolute_error /
        std::max({std::abs(central), std::abs(analytic), 1e-24L});
    return result;
}

void write_value(
    const char* name,
    const LocalSldProjectiveHeightGapCorrelationObjectiveValue& value,
    std::ostream& output,
    bool trailing_comma) {
    output << std::setprecision(18)
        << "  \"" << name << "\": {"
        << "\"gram_pairing\": "
        << static_cast<double>(value.gram_pairing)
        << ", \"first_h2_norm2\": "
        << static_cast<double>(value.first_h2_norm2)
        << ", \"second_h2_norm2\": "
        << static_cast<double>(value.second_h2_norm2)
        << ", \"correlation_squared\": "
        << static_cast<double>(value.correlation_squared)
        << ", \"half_decay_weighted_correlation\": "
        << static_cast<double>(value.half_decay_weighted_correlation)
        << ", \"weighted_correlation_squared\": "
        << static_cast<double>(value.weighted_correlation_squared)
        << ", \"finite\": " << (value.finite ? "true" : "false")
        << '}' << (trailing_comma ? ",\n" : "\n");
}

}  // namespace

LocalSldProjectiveHeightGapCorrelationAdversaryOptions
LocalSldProjectiveHeightGapCorrelationAdversaryCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveHeightGapCorrelationAdversaryOptions options;
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
        } else if (name == "--first-shell") {
            options.first_shell = std::stoi(next(index, name));
        } else if (name == "--second-shell") {
            options.second_shell = std::stoi(next(index, name));
        } else if (name == "--iterations") {
            options.iterations = std::stoi(next(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next(index, name));
        } else if (name == "--lbfgs-history") {
            options.lbfgs_history = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next(index, name));
        } else if (name == "--validation-seed") {
            options.validation_seed = std::stoull(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown height-gap-correlation option: " + name);
        }
    }
    if (options.state_path.empty() ||
        options.output_state_path.empty() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.first_shell < 0 ||
        options.second_shell <= options.first_shell ||
        options.second_shell > 30 ||
        options.iterations < 0 || options.iterations > 10000 ||
        options.line_search_steps < 1 ||
        options.lbfgs_history < 1 || options.lbfgs_history > 64 ||
        options.threads < 1 || options.threads > 256 ||
        !(options.initial_step > 0.0L) ||
        !std::isfinite(options.initial_step)) {
        throw std::invalid_argument(
            "height-gap-correlation requires state/output/certificate, valid selection, 0<=first<second<=30, nonnegative iterations, and valid search parameters");
    }
    return options;
}

void LocalSldProjectiveHeightGapCorrelationAdversaryCli::print_help(
    std::ostream& out) {
    out << "Projective height-gap correlation adversary options:\n"
        << "  --state PATH           initial Fourier TSV\n"
        << "  --output-state PATH    optimized Fourier TSV\n"
        << "  --certificate PATH     write English JSON trace\n"
        << "  --selection NAME       local SLD triad selection\n"
        << "  --first-shell I        lower primitive-height shell\n"
        << "  --second-shell J       upper primitive-height shell\n"
        << "  --iterations N         exact-gradient L-BFGS steps\n"
        << "  --line-search N        trials per step\n"
        << "  --lbfgs-history N      curvature pairs\n"
        << "  --threads N            direct workers\n"
        << "  --step X               initial line-search step\n"
        << "  --validation-seed N    deterministic gradient check\n";
}

int LocalSldProjectiveHeightGapCorrelationAdversaryCli::run(
    const LocalSldProjectiveHeightGapCorrelationAdversaryOptions& options,
    std::ostream& out) {
    const SpectralState initial = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const SpectralObjective spectral_objective(dynamics);
    const SpectralAdjoint adjoint(dynamics, spectral_objective);
    const GradientAdversary adversary(
        dynamics, spectral_objective, adjoint);
    const TriadSelection selection = LocalSldTriadSelection::parse(
        options.selection);
    const LocalSldProjectiveHeightGapCorrelationObjective objective(
        selection, options.first_shell, options.second_shell,
        options.threads);
    const auto initial_value = objective.evaluate(initial);
    const GradientDirectionalCheck initial_gradient =
        gradient_directional_check(
        dynamics, objective, initial, options.validation_seed);

    GradientSearchOptions search;
    search.iterations = options.iterations;
    search.line_search_steps = options.line_search_steps;
    search.initial_step = options.initial_step;
    search.method = "lbfgs";
    search.lbfgs_history = options.lbfgs_history;
    search.objective =
        "local-projective-height-gap-correlation-ratio";
    search.closure_selection = selection;
    search.projective_first_height_shell = options.first_shell;
    search.projective_second_height_shell = options.second_shell;
    search.objective_threads = options.threads;
    const GradientSearchResult optimized = adversary.maximize_q(
        initial, search);
    const auto final_value = objective.evaluate(optimized.state);
    const GradientDirectionalCheck final_gradient =
        gradient_directional_check(
            dynamics, objective, optimized.state,
            options.validation_seed + 1U);

    std::ostringstream metadata;
    metadata << std::setprecision(18)
        << "projective height-gap correlation adversary; shells="
        << options.first_shell << ',' << options.second_shell
        << "; weighted_correlation="
        << static_cast<double>(
               final_value.half_decay_weighted_correlation)
        << "; candidate_lemma_proved=false";
    const std::filesystem::path state_path(options.output_state_path);
    if (!state_path.parent_path().empty()) {
        std::filesystem::create_directories(state_path.parent_path());
    }
    SpectralStateWriter::write_tsv(
        options.output_state_path, optimized.state, metadata.str());

    const std::filesystem::path certificate_path(options.certificate_path);
    if (!certificate_path.parent_path().empty()) {
        std::filesystem::create_directories(
            certificate_path.parent_path());
    }
    std::ofstream certificate(certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write height-gap-correlation certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-height-gap-correlation-adversary-v1\",\n"
        << "  \"initial_state_path\": \"" << options.state_path
        << "\",\n"
        << "  \"state_path\": \"" << options.output_state_path
        << "\",\n"
        << "  \"triad_selection\": \"" << options.selection
        << "\",\n"
        << "  \"first_shell\": " << options.first_shell << ",\n"
        << "  \"second_shell\": " << options.second_shell << ",\n"
        << "  \"shell_gap\": "
        << options.second_shell - options.first_shell << ",\n"
        << "  \"iterations_requested\": " << options.iterations << ",\n"
        << "  \"iterations_executed\": " << optimized.iterations << ",\n"
        << "  \"accepted_steps\": " << optimized.accepted_steps << ",\n"
        << "  \"evaluations\": " << optimized.trajectory_evaluations
        << ",\n"
        << "  \"final_projected_gradient_norm\": "
        << static_cast<double>(optimized.final_projected_gradient_norm)
        << ",\n"
        << "  \"gradient_directional_relative_error\": "
        << static_cast<double>(initial_gradient.relative_error) << ",\n"
        << "  \"initial_gradient_directional_absolute_error\": "
        << static_cast<double>(initial_gradient.absolute_error) << ",\n"
        << "  \"initial_gradient_directional_relative_error\": "
        << static_cast<double>(initial_gradient.relative_error) << ",\n"
        << "  \"final_gradient_directional_absolute_error\": "
        << static_cast<double>(final_gradient.absolute_error) << ",\n"
        << "  \"final_gradient_directional_relative_error\": "
        << static_cast<double>(final_gradient.relative_error) << ",\n";
    write_value("initial", initial_value, certificate, true);
    write_value("final", final_value, certificate, true);
    certificate
        << "  \"finite_search_is_not_a_proof\": true,\n"
        << "  \"uniform_half_gap_decay_proved\": false,\n"
        << "  \"interpretation\": \"exact-gradient falsification search for the PNT-13 normalized primitive-height gap estimate\"\n"
        << "}\n";
    out << std::setprecision(12)
        << "projective height-gap correlation shells="
        << options.first_shell << ',' << options.second_shell
        << " weighted="
        << static_cast<double>(
               final_value.half_decay_weighted_correlation)
        << " gain="
        << static_cast<double>(
               final_value.weighted_correlation_squared /
               std::max(initial_value.weighted_correlation_squared, 1e-30L))
        << " initial_gradient_error="
        << static_cast<double>(initial_gradient.relative_error)
        << " final_gradient_error="
        << static_cast<double>(final_gradient.relative_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return initial_value.finite && final_value.finite &&
            std::isfinite(initial_gradient.relative_error) &&
            std::isfinite(final_gradient.relative_error)
        ? 0 : 2;
}

}  // namespace lemma
