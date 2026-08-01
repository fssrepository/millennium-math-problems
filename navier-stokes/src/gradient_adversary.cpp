#include "gradient_adversary.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lemma {
namespace {

SpectralReal increment_inner_product(const SpectralIncrement& left,
                                     const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("gradient inner-product layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal project_to_energy_sphere(SpectralIncrement& gradient,
                                      const SpectralState& state) {
    const SpectralReal energy = SpectralStateOps::energy(state);
    if (!(energy > 0.0L)) {
        throw std::invalid_argument("cannot optimize a zero-energy state");
    }
    const SpectralReal radial_coefficient =
        increment_inner_product(gradient, state.velocity) / energy;
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            gradient[mode][component] -=
                radial_coefficient * state.velocity[mode][component];
        }
    }
    const SpectralReal norm2 = increment_inner_product(gradient, gradient);
    return std::sqrt(std::max(0.0L, norm2));
}

SpectralReal increment_norm(const SpectralIncrement& increment) {
    return std::sqrt(std::max(
        0.0L, increment_inner_product(increment, increment)));
}

}  // namespace

GradientAdversary::GradientAdversary(const SpectralDynamics& dynamics,
                                     const SpectralObjective& objective,
                                     const SpectralAdjoint& adjoint)
    : dynamics_(dynamics), objective_(objective), adjoint_(adjoint) {}

SpectralReal GradientAdversary::objective_value(
    const SpectralState& initial,
    const GradientSearchOptions& options) const {
    SpectralState state = initial;
    const SpectralReal initial_q =
        objective_.evaluate(state).energy_level_quantity;
    SpectralReal previous_integrand =
        objective_.evaluate(state).critical_integrand;
    SpectralReal critical_integral = 0.0L;
    SpectralReal maximum_q = initial_q;
    for (int step = 0; step < options.trajectory_steps; ++step) {
        dynamics_.rk4_step(state, options.viscosity, options.time_step);
        const StaticObjective sample = objective_.evaluate(state);
        maximum_q = std::max(maximum_q, sample.energy_level_quantity);
        critical_integral += 0.5L * options.time_step *
                             (previous_integrand + sample.critical_integrand);
        previous_integrand = sample.critical_integrand;
    }
    const SpectralReal terminal_q =
        objective_.evaluate(state).energy_level_quantity;
    if (options.objective == "max-q") {
        return maximum_q;
    }
    if (options.objective == "terminal-q") {
        return terminal_q;
    }
    if (options.objective == "q-gain") {
        if (!(initial_q > 1e-30L) || !(terminal_q > 1e-30L)) {
            return -std::numeric_limits<SpectralReal>::infinity();
        }
        return std::log(terminal_q / initial_q);
    }
    if (options.objective == "q-increase") {
        return terminal_q - initial_q;
    }
    if (options.objective == "critical-integral") {
        return critical_integral;
    }
    throw std::invalid_argument("unknown gradient objective: " +
                                options.objective);
}

GradientSearchResult GradientAdversary::maximize_q(
    const SpectralState& initial,
    const GradientSearchOptions& options) const {
    if (options.iterations < 0 || options.line_search_steps < 1 ||
        options.trajectory_steps < 0) {
        throw std::invalid_argument("invalid gradient-search iteration count");
    }
    if (!(options.initial_step > 0.0L)) {
        throw std::invalid_argument("gradient-search step must be positive");
    }
    GradientSearchResult result;
    const InitialSobolevConstraint sobolev(
        options.sobolev_order, options.sobolev_cap);
    result.state = initial;
    const SpectralReal target_energy = SpectralStateOps::energy(initial);
    dynamics_.enforce_constraints(result.state);
    SpectralStateOps::normalize_energy(result.state, target_energy);
    if (!sobolev.admissible(result.state)) {
        throw std::invalid_argument(
            "initial state exceeds the configured Sobolev cap");
    }
    result.objective = objective_value(result.state, options);
    result.initial_objective = result.objective;
    ++result.trajectory_evaluations;
    SpectralReal next_step = options.initial_step;

    for (int iteration = 0; iteration < options.iterations; ++iteration) {
        QTrajectoryGradient trajectory;
        if (options.objective == "max-q") {
            trajectory = adjoint_.maximum_q_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps);
        } else if (options.objective == "terminal-q") {
            trajectory = adjoint_.terminal_q_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps);
        } else if (options.objective == "q-gain") {
            trajectory = adjoint_.q_gain_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps);
        } else if (options.objective == "q-increase") {
            trajectory = adjoint_.q_increase_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps);
        } else if (options.objective == "critical-integral") {
            trajectory = adjoint_.critical_integral_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps);
        } else {
            throw std::invalid_argument("unknown gradient objective: " +
                                        options.objective);
        }
        ++result.trajectory_evaluations;
        ++result.iterations;
        result.objective = trajectory.objective_value;
        result.objective_step = trajectory.objective_step;
        GradientIterationRecord record;
        record.iteration = iteration;
        record.objective_before = result.objective;
        SpectralIncrement direction = trajectory.initial_gradient;
        SpectralReal gradient_norm =
            project_to_energy_sphere(direction, result.state);
        const SpectralReal sobolev_value = sobolev.value(result.state);
        if (sobolev.enabled() &&
            sobolev_value >= 0.99L * sobolev.cap()) {
            const SpectralIncrement normal =
                sobolev.energy_tangent_normal(result.state);
            const SpectralReal normal_norm2 =
                increment_inner_product(normal, normal);
            const SpectralReal outward =
                increment_inner_product(direction, normal);
            if (outward > 0.0L && normal_norm2 > 1e-30L) {
                const SpectralReal coefficient = outward / normal_norm2;
                for (std::size_t mode = 0; mode < direction.size(); ++mode) {
                    for (std::size_t component = 0; component < 3;
                         ++component) {
                        direction[mode][component] -=
                            coefficient * normal[mode][component];
                    }
                }
                gradient_norm = increment_norm(direction);
            }
        }
        result.final_projected_gradient_norm = gradient_norm;
        record.projected_gradient_norm = gradient_norm;
        record.sobolev_value = sobolev_value;
        if (!(gradient_norm > 1e-24L) || !std::isfinite(gradient_norm)) {
            record.objective_after = result.objective;
            result.trace.push_back(record);
            break;
        }
        for (ComplexVector& value : direction) {
            for (SpectralComplex& component : value) {
                component /= gradient_norm;
            }
        }

        bool accepted = false;
        SpectralReal trial_step = next_step;
        for (int line_step = 0;
             line_step < options.line_search_steps; ++line_step) {
            ++record.line_search_evaluations;
            SpectralState candidate = dynamics_.add_increment(
                result.state, direction, trial_step);
            dynamics_.enforce_constraints(candidate);
            sobolev.retract(candidate, target_energy);
            const SpectralReal candidate_objective =
                objective_value(candidate, options);
            ++result.trajectory_evaluations;
            const SpectralReal improvement_floor =
                1e-13L * std::max(1e-30L, std::abs(result.objective));
            if (std::isfinite(candidate_objective) &&
                candidate_objective > result.objective + improvement_floor) {
                result.state = std::move(candidate);
                result.objective = candidate_objective;
                ++result.accepted_steps;
                next_step = std::min(
                    options.initial_step, 1.5L * trial_step);
                record.accepted_step = trial_step;
                record.accepted = true;
                accepted = true;
                break;
            }
            trial_step *= 0.5L;
        }
        record.objective_after = result.objective;
        result.trace.push_back(record);
        if (!accepted) {
            break;
        }
    }
    result.final_sobolev_value = sobolev.value(result.state);
    return result;
}

}  // namespace lemma
