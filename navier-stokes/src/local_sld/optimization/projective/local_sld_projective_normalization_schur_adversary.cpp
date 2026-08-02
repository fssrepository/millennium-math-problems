#include "local_sld_projective_normalization_schur_adversary.hpp"

#include "gradient_adversary.hpp"
#include "local_sld_projective_normalization_schur_objective.hpp"
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

struct GradientCheck {
    SpectralReal absolute_error = 0.0L;
    SpectralReal relative_error = 0.0L;
};

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

GradientCheck gradient_check(
    const SpectralDynamics& dynamics,
    const LocalSldProjectiveNormalizationSchurObjective& objective,
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
    SpectralState plus = dynamics.add_increment(state, direction, epsilon);
    SpectralState minus = dynamics.add_increment(state, direction, -epsilon);
    dynamics.enforce_constraints(plus);
    dynamics.enforce_constraints(minus);
    const SpectralReal central =
        (objective.evaluate(plus)
             .height_half_compensated_schur_squared_majorant -
         objective.evaluate(minus)
             .height_half_compensated_schur_squared_majorant) /
        (2.0L * epsilon);
    GradientCheck result;
    result.absolute_error = std::abs(central - analytic);
    result.relative_error = result.absolute_error /
        std::max({std::abs(central), std::abs(analytic), 1e-30L});
    return result;
}

void write_value(
    const char* name,
    const LocalSldProjectiveNormalizationSchurObjectiveValue& value,
    std::ostream& output,
    bool trailing_comma) {
    output << std::setprecision(18)
        << "  \"" << name << "\": {"
        << "\"full_stretching\": "
        << static_cast<double>(value.full_stretching)
        << ", \"enstrophy\": "
        << static_cast<double>(value.enstrophy)
        << ", \"palinstrophy\": "
        << static_cast<double>(value.palinstrophy)
        << ", \"selected_aggregate_h1_norm2\": "
        << static_cast<double>(value.selected_aggregate_h1_norm2)
        << ", \"diagonal_tail_h2_norm2\": "
        << static_cast<double>(value.diagonal_tail_h2_norm2)
        << ", \"normalized_absolute_gram_row_sum\": "
        << static_cast<double>(value.normalized_absolute_gram_row_sum)
        << ", \"normalization_common_factor\": "
        << static_cast<double>(value.normalization_common_factor)
        << ", \"schur_squared_majorant\": "
        << static_cast<double>(value.schur_squared_majorant)
        << ", \"height_half_compensated_schur_squared_majorant\": "
        << static_cast<double>(
               value.height_half_compensated_schur_squared_majorant)
        << ", \"finite\": " << (value.finite ? "true" : "false")
        << '}' << (trailing_comma ? ",\n" : "\n");
}

}  // namespace

LocalSldProjectiveNormalizationSchurAdversaryOptions
LocalSldProjectiveNormalizationSchurAdversaryCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveNormalizationSchurAdversaryOptions options;
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
        } else if (name == "--projective-core-height") {
            options.core_maximum_height = std::stoll(next(index, name));
        } else if (name == "--row-shell") {
            options.row_shell = std::stoi(next(index, name));
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
                "unknown normalization-Schur option: " + name);
        }
    }
    if (options.state_path.empty() || options.output_state_path.empty() ||
        options.certificate_path.empty() ||
        !LocalSldTriadSelection::supports(options.selection) ||
        options.core_maximum_height < 1 || options.row_shell < 0 ||
        options.row_shell > 62 || options.iterations < 0 ||
        options.iterations > 10000 || options.line_search_steps < 1 ||
        options.lbfgs_history < 1 || options.lbfgs_history > 64 ||
        options.threads < 1 || options.threads > 256 ||
        !(options.initial_step > 0.0L) ||
        !std::isfinite(options.initial_step)) {
        throw std::invalid_argument(
            "normalization-Schur requires state/output/certificate, valid selection/height/row, nonnegative iterations, and valid search parameters");
    }
    return options;
}

void LocalSldProjectiveNormalizationSchurAdversaryCli::print_help(
    std::ostream& out) {
    out << "Projective normalization Schur adversary options:\n"
        << "  --state PATH                 initial Fourier TSV\n"
        << "  --output-state PATH          optimized Fourier TSV\n"
        << "  --certificate PATH           write English JSON trace\n"
        << "  --selection NAME             local SLD triad selection\n"
        << "  --projective-core-height H   tail threshold\n"
        << "  --row-shell I                differentiated Schur row\n"
        << "  --iterations N               exact-gradient L-BFGS steps\n"
        << "  --line-search N              trials per step\n"
        << "  --lbfgs-history N            curvature pairs\n"
        << "  --threads N                  direct workers\n"
        << "  --step X                     initial line-search step\n"
        << "  --validation-seed N          deterministic gradient check\n";
}

int LocalSldProjectiveNormalizationSchurAdversaryCli::run(
    const LocalSldProjectiveNormalizationSchurAdversaryOptions& options,
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
    const LocalSldProjectiveNormalizationSchurObjective objective(
        dynamics, selection, options.core_maximum_height,
        options.row_shell, options.threads);
    const auto initial_value = objective.evaluate(initial);
    const GradientCheck initial_gradient = gradient_check(
        dynamics, objective, initial, options.validation_seed);

    GradientSearchOptions search;
    search.iterations = options.iterations;
    search.line_search_steps = options.line_search_steps;
    search.initial_step = options.initial_step;
    search.method = "lbfgs";
    search.lbfgs_history = options.lbfgs_history;
    search.objective =
        "local-projective-normalization-schur-row-ratio";
    search.closure_selection = selection;
    search.projective_core_maximum_height =
        options.core_maximum_height;
    search.projective_schur_row_shell = options.row_shell;
    search.objective_threads = options.threads;
    const GradientSearchResult optimized = adversary.maximize_q(
        initial, search);
    const auto final_value = objective.evaluate(optimized.state);
    const GradientCheck final_gradient = gradient_check(
        dynamics, objective, optimized.state,
        options.validation_seed + 1U);

    std::ostringstream metadata;
    metadata << std::setprecision(18)
        << "fixed-row PNT-12 Schur adversary; H="
        << options.core_maximum_height << "; row=" << options.row_shell
        << "; compensated_schur="
        << static_cast<double>(
               final_value.height_half_compensated_schur_squared_majorant)
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
            "cannot write normalization-Schur certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-normalization-schur-adversary-v1\",\n"
        << "  \"initial_state_path\": \"" << options.state_path
        << "\",\n"
        << "  \"state_path\": \"" << options.output_state_path
        << "\",\n"
        << "  \"triad_selection\": \"" << options.selection
        << "\",\n"
        << "  \"core_maximum_height\": "
        << options.core_maximum_height << ",\n"
        << "  \"row_shell\": " << options.row_shell << ",\n"
        << "  \"iterations_requested\": " << options.iterations << ",\n"
        << "  \"iterations_executed\": " << optimized.iterations << ",\n"
        << "  \"accepted_steps\": " << optimized.accepted_steps << ",\n"
        << "  \"evaluations\": " << optimized.trajectory_evaluations
        << ",\n"
        << "  \"final_projected_gradient_norm\": "
        << static_cast<double>(optimized.final_projected_gradient_norm)
        << ",\n"
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
        << "  \"trace\": [\n";
    for (std::size_t index = 0; index < optimized.trace.size(); ++index) {
        const auto& row = optimized.trace[index];
        certificate
            << "    {\"iteration\": " << row.iteration
            << ", \"objective_before\": "
            << static_cast<double>(row.objective_before)
            << ", \"objective_after\": "
            << static_cast<double>(row.objective_after)
            << ", \"projected_gradient_norm\": "
            << static_cast<double>(row.projected_gradient_norm)
            << ", \"accepted_step\": "
            << static_cast<double>(row.accepted_step)
            << ", \"line_search_evaluations\": "
            << row.line_search_evaluations
            << ", \"used_steepest_fallback\": "
            << (row.used_steepest_fallback ? "true" : "false")
            << ", \"accepted\": "
            << (row.accepted ? "true" : "false") << '}'
            << (index + 1 == optimized.trace.size() ? "\n" : ",\n");
    }
    certificate
        << "  ],\n"
        << "  \"finite_search_is_not_a_proof\": true,\n"
        << "  \"uniform_PNT12_bound_proved\": false,\n"
        << "  \"interpretation\": \"exact-gradient stress search for one smooth fixed-row branch of the coupled PNT-12 Schur majorant\"\n"
        << "}\n";
    out << std::setprecision(12)
        << "normalization Schur H=" << options.core_maximum_height
        << " row=" << options.row_shell
        << " objective=" << static_cast<double>(
               final_value.height_half_compensated_schur_squared_majorant)
        << " gain=" << static_cast<double>(
               final_value.height_half_compensated_schur_squared_majorant /
               std::max(
                   initial_value
                       .height_half_compensated_schur_squared_majorant,
                   1e-30L))
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
