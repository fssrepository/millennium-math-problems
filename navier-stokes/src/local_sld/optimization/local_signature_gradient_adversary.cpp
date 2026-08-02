#include "local_signature_gradient_adversary.hpp"

#include "local_signature_objective.hpp"
#include "parallel_executor.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

SpectralReal inner_product(const SpectralIncrement& left,
                           const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "local signature gradient layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

void project_spectral_constraints(
    SpectralIncrement& gradient, const SpectralState& state) {
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const WaveVector wave = state.waves[mode];
        if (!is_positive_representative(wave)) {
            continue;
        }
        const std::size_t negative = state.index.at(-wave);
        ComplexVector combined{};
        const ComplexVector reflected = conjugate(gradient[negative]);
        for (std::size_t component = 0; component < 3; ++component) {
            combined[component] = 0.5L *
                (gradient[mode][component] + reflected[component]);
        }
        combined = project_divergence_free(wave, combined);
        gradient[mode] = combined;
        gradient[negative] = conjugate(combined);
    }
}

SpectralReal project_to_energy_sphere(
    SpectralIncrement& gradient, const SpectralState& state) {
    project_spectral_constraints(gradient, state);
    const SpectralReal energy = SpectralStateOps::energy(state);
    if (!(energy > 0.0L)) {
        throw std::invalid_argument("cannot optimize zero-energy state");
    }
    const SpectralReal radial =
        inner_product(gradient, state.velocity) / energy;
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            gradient[mode][component] -=
                radial * state.velocity[mode][component];
        }
    }
    return std::sqrt(std::max(0.0L, inner_product(gradient, gradient)));
}

SpectralState retract(
    const SpectralState& state, const SpectralIncrement& direction,
    SpectralReal step, SpectralReal energy) {
    SpectralState candidate = state;
    for (std::size_t mode = 0; mode < candidate.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            candidate.velocity[mode][component] +=
                step * direction[mode][component];
        }
    }
    SpectralStateOps::normalize_energy(candidate, energy);
    return candidate;
}

SpectralReal objective_value(
    const SpectralState& state, const std::string& objective) {
    const LocalSignatureObjectiveValue value =
        LocalSignatureObjective::evaluate(state);
    if (objective == "amplification") {
        return value.signed_amplification;
    }
    if (objective == "transfer") {
        return std::abs(value.signed_local_transfer);
    }
    throw std::invalid_argument(
        "local signature objective must be amplification or transfer");
}

SpectralIncrement objective_gradient(
    const SpectralState& state, const std::string& objective) {
    if (objective == "amplification") {
        return LocalSignatureObjective::signed_amplification_gradient(state);
    }
    if (objective == "transfer") {
        return LocalSignatureObjective::absolute_signed_transfer_gradient(
            state);
    }
    throw std::invalid_argument(
        "local signature objective must be amplification or transfer");
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

SpectralReal fitted_log_slope(
    const std::vector<LocalSignatureGradientCutoffRow>& rows) {
    SpectralReal count = 0.0L;
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    for (const auto& row : rows) {
        if (row.cutoff <= 1 || !(row.best_objective > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
        const SpectralReal y = std::log(row.best_objective);
        count += 1.0L;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }
    const SpectralReal denominator = count * sum_xx - sum_x * sum_x;
    if (count < 2.0L || std::abs(denominator) <
            std::numeric_limits<SpectralReal>::epsilon()) {
        return 0.0L;
    }
    return (count * sum_xy - sum_x * sum_y) / denominator;
}

}  // namespace

LocalSignatureGradientResult LocalSignatureGradientAdversary::maximize(
    const SpectralState& initial,
    const LocalSignatureGradientOptions& options) {
    if (options.iterations < 1 || options.line_search_steps < 1 ||
        !(options.initial_step > 0.0L) ||
        (options.objective != "amplification" &&
         options.objective != "transfer")) {
        throw std::invalid_argument(
            "invalid local signature gradient adversary options");
    }
    LocalSignatureGradientResult result;
    result.state = initial;
    const SpectralReal energy = SpectralStateOps::energy(initial);
    result.objective = objective_value(result.state, options.objective);
    result.initial_objective = result.objective;
    ++result.evaluations;
    for (int iteration = 0; iteration < options.iterations; ++iteration) {
        SpectralIncrement direction = objective_gradient(
            result.state, options.objective);
        const SpectralReal gradient_norm = project_to_energy_sphere(
            direction, result.state);
        LocalSignatureGradientTraceRow row;
        row.iteration = iteration;
        row.objective = result.objective;
        row.projected_gradient_norm = gradient_norm;
        if (!(gradient_norm > 0.0L) || !std::isfinite(gradient_norm)) {
            result.trace.push_back(row);
            break;
        }
        for (ComplexVector& value : direction) {
            for (SpectralComplex& component : value) {
                component /= gradient_norm;
            }
        }
        SpectralReal step = options.initial_step;
        for (int line = 0; line < options.line_search_steps; ++line) {
            SpectralState candidate = retract(
                result.state, direction, step, energy);
            const SpectralReal candidate_objective = objective_value(
                candidate, options.objective);
            ++result.evaluations;
            if (candidate_objective > result.objective) {
                result.state = std::move(candidate);
                result.objective = candidate_objective;
                row.accepted_step = step;
                row.accepted = true;
                ++result.accepted_steps;
                break;
            }
            step *= 0.5L;
        }
        result.trace.push_back(row);
        if (!row.accepted) {
            break;
        }
    }
    return result;
}

LocalSignatureGradientCliOptions LocalSignatureGradientCli::parse(
    int argc, char** argv, int first) {
    LocalSignatureGradientCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--min-cutoff") {
            options.minimum_cutoff = std::stoi(next(index, name));
        } else if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--restarts") {
            options.restarts = std::stoi(next(index, name));
        } else if (name == "--workers") {
            options.workers = std::stoi(next(index, name));
        } else if (name == "--iterations") {
            options.iterations = std::stoi(next(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next(index, name));
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else if (name == "--profile") {
            options.profile = next(index, name);
        } else if (name == "--objective") {
            options.objective = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-signature-gradient option: " + name);
        }
    }
    return options;
}

void LocalSignatureGradientCli::print_help(std::ostream& out) {
    out << "Local-signature gradient adversary options:\n"
        << "  --min-cutoff K       first Fourier cutoff\n"
        << "  --max-cutoff K       last Fourier cutoff\n"
        << "  --profile NAME       decaying, flat, or outer-half-flat\n"
        << "  --objective NAME     amplification or transfer\n"
        << "  --restarts N         independent exact-gradient starts\n"
        << "  --workers N          parallel restart workers\n"
        << "  --iterations N       gradient iterations per restart\n"
        << "  --line-search N      backtracking trials per iteration\n"
        << "  --step X             initial Riemannian step\n"
        << "  --seed N             deterministic master seed\n"
        << "  --certificate PATH   write optimization JSON\n";
}

int LocalSignatureGradientCli::run(
    const LocalSignatureGradientCliOptions& options, std::ostream& out) {
    if (options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 8 || options.restarts < 1 ||
        options.workers < 1 || options.iterations < 1 ||
        options.certificate_path.empty() ||
        (options.objective != "amplification" &&
         options.objective != "transfer")) {
        throw std::invalid_argument(
            "local-signature-gradient requires valid cutoffs, positive search counts, and --certificate");
    }
    const LocalSignatureStateProfile profile =
        LocalSignatureStateFactory::parse(options.profile);
    const std::size_t cutoff_count = static_cast<std::size_t>(
        options.maximum_cutoff - options.minimum_cutoff + 1);
    const std::size_t task_count = cutoff_count *
        static_cast<std::size_t>(options.restarts);
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        const SpectralState warm = LocalSignatureStateFactory::make(
            cutoff, profile, options.seed);
        static_cast<void>(SpectralStateOps::interactions(warm));
    }
    std::vector<LocalSignatureGradientResult> results(task_count);
    std::vector<std::uint64_t> seeds(task_count);
    LocalSignatureGradientOptions search;
    search.iterations = options.iterations;
    search.objective = options.objective;
    search.line_search_steps = options.line_search_steps;
    search.initial_step = options.initial_step;
    const ParallelExecutor executor(options.workers);
    executor.for_each(task_count, [&](std::size_t task) {
        const int cutoff = options.minimum_cutoff + static_cast<int>(
            task / static_cast<std::size_t>(options.restarts));
        const std::size_t restart = task %
            static_cast<std::size_t>(options.restarts);
        const std::uint64_t seed = splitmix64(
            options.seed ^ splitmix64(static_cast<std::uint64_t>(cutoff)) ^
            splitmix64(static_cast<std::uint64_t>(restart + 1)));
        seeds[task] = seed;
        results[task] = LocalSignatureGradientAdversary::maximize(
            LocalSignatureStateFactory::make(cutoff, profile, seed), search);
    });

    LocalSignatureGradientReport report;
    report.profile = LocalSignatureStateFactory::name(profile);
    report.objective = options.objective;
    report.workers = executor.threads();
    report.restarts = options.restarts;
    report.iterations = options.iterations;
    report.critical_log_slope = options.objective == "amplification"
        ? 0.5L
        : 4.0L;
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        LocalSignatureGradientCutoffRow row;
        row.cutoff = cutoff;
        const std::size_t base = static_cast<std::size_t>(
            cutoff - options.minimum_cutoff) *
            static_cast<std::size_t>(options.restarts);
        for (int restart = 0; restart < options.restarts; ++restart) {
            const std::size_t index = base +
                static_cast<std::size_t>(restart);
            const auto& candidate = results[index];
            if (candidate.objective > row.best_objective) {
                const LocalSignatureObjectiveValue value =
                    LocalSignatureObjective::evaluate(candidate.state);
                const SpectralReal frequency =
                    static_cast<SpectralReal>(cutoff);
                const SpectralReal energy =
                    SpectralStateOps::energy(candidate.state);
                row.best_initial_objective = candidate.initial_objective;
                row.best_objective = candidate.objective;
                row.absolute_signed_transfer =
                    std::abs(value.signed_local_transfer);
                row.signature_transfer_l2 =
                    std::sqrt(value.squared_signature_transfer);
                row.normalized_lsf2_l2 = row.signature_transfer_l2 /
                    (std::pow(frequency, 3.5L) *
                     std::pow(energy, 1.5L));
                row.normalized_viscous_transfer =
                    row.absolute_signed_transfer /
                    (std::pow(frequency, 4.0L) * energy);
                row.best_seed = seeds[index];
                row.accepted_steps = candidate.accepted_steps;
                row.evaluations = candidate.evaluations;
            }
        }
        if (row.best_initial_objective > 0.0L) {
            row.improvement_factor =
                row.best_objective / row.best_initial_objective;
        }
        report.rows.push_back(row);
    }
    report.fitted_log_slope = fitted_log_slope(report.rows);

    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write local signature gradient certificate");
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-local-signature-gradient-v1\",\n"
        << "  \"profile\": \"" << report.profile
        << "\",\n  \"objective\": \"" << report.objective
        << "\",\n  \"workers\": " << report.workers
        << ",\n  \"restarts\": " << report.restarts
        << ",\n  \"iterations\": " << report.iterations
        << ",\n  \"fitted_log_slope\": "
        << static_cast<double>(report.fitted_log_slope)
        << ",\n  \"critical_log_slope\": "
        << static_cast<double>(report.critical_log_slope) << ",\n"
        << "  \"finite_search_is_not_a_proof\": true,\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"best_initial_objective\": "
            << static_cast<double>(row.best_initial_objective)
            << ", \"best_objective\": "
            << static_cast<double>(row.best_objective)
            << ", \"improvement_factor\": "
            << static_cast<double>(row.improvement_factor)
            << ", \"absolute_signed_transfer\": "
            << static_cast<double>(row.absolute_signed_transfer)
            << ", \"signature_transfer_l2\": "
            << static_cast<double>(row.signature_transfer_l2)
            << ", \"normalized_lsf2_l2\": "
            << static_cast<double>(row.normalized_lsf2_l2)
            << ", \"normalized_viscous_transfer\": "
            << static_cast<double>(row.normalized_viscous_transfer)
            << ", \"best_seed\": \"" << row.best_seed
            << "\", \"accepted_steps\": " << row.accepted_steps
            << ", \"evaluations\": " << row.evaluations << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << "local signature exact-gradient adversary profile="
        << report.profile << " objective=" << report.objective
        << " workers=" << report.workers
        << " restarts=" << report.restarts
        << "\ncutoff,initial_objective,optimized_objective,improvement,"
           "signed_transfer,signature_l2,lsf2_ratio,viscous_ratio,seed,"
           "accepted,evaluations\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ','
            << static_cast<double>(row.best_initial_objective) << ','
            << static_cast<double>(row.best_objective) << ','
            << static_cast<double>(row.improvement_factor) << ','
            << static_cast<double>(row.absolute_signed_transfer) << ','
            << static_cast<double>(row.signature_transfer_l2) << ','
            << static_cast<double>(row.normalized_lsf2_l2) << ','
            << static_cast<double>(row.normalized_viscous_transfer) << ','
            << row.best_seed << ',' << row.accepted_steps << ','
            << row.evaluations << '\n';
    }
    out << "fitted optimized slope="
        << static_cast<double>(report.fitted_log_slope)
        << " critical=" << static_cast<double>(report.critical_log_slope)
        << "\nCertificate written to "
        << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
