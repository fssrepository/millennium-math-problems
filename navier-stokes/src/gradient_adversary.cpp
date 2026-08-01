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

void add_scaled_increment(SpectralIncrement& target,
                          const SpectralIncrement& source,
                          SpectralReal scale) {
    if (target.size() != source.size()) {
        throw std::invalid_argument("L-BFGS increment layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += scale * source[mode][component];
        }
    }
}

SpectralIncrement increment_difference(const SpectralIncrement& left,
                                       const SpectralIncrement& right) {
    SpectralIncrement result = left;
    add_scaled_increment(result, right, -1.0L);
    return result;
}

SpectralReal project_to_search_tangent(
    SpectralIncrement& increment, const SpectralState& state,
    const InitialSobolevConstraint& sobolev) {
    project_to_energy_sphere(increment, state);
    const SpectralReal sobolev_value = sobolev.value(state);
    if (sobolev.enabled() && sobolev_value >= 0.99L * sobolev.cap()) {
        const SpectralIncrement normal =
            sobolev.energy_tangent_normal(state);
        const SpectralReal normal_norm2 =
            increment_inner_product(normal, normal);
        const SpectralReal outward =
            increment_inner_product(increment, normal);
        if (outward > 0.0L && normal_norm2 > 1e-30L) {
            add_scaled_increment(
                increment, normal, -outward / normal_norm2);
        }
    }
    return increment_norm(increment);
}

struct LbfgsPair {
    SpectralIncrement state_delta;
    SpectralIncrement gradient_delta;
    SpectralReal inverse_curvature = 0.0L;
};

SpectralIncrement lbfgs_ascent_direction(
    const SpectralIncrement& gradient,
    const std::vector<LbfgsPair>& history) {
    if (history.empty()) {
        return gradient;
    }
    SpectralIncrement direction = gradient;
    std::vector<SpectralReal> alpha(history.size(), 0.0L);
    for (std::size_t reverse = history.size(); reverse > 0; --reverse) {
        const std::size_t index = reverse - 1;
        const LbfgsPair& pair = history[index];
        alpha[index] = pair.inverse_curvature *
                       increment_inner_product(pair.state_delta, direction);
        add_scaled_increment(
            direction, pair.gradient_delta, -alpha[index]);
    }
    const LbfgsPair& newest = history.back();
    const SpectralReal yy = increment_inner_product(
        newest.gradient_delta, newest.gradient_delta);
    const SpectralReal sy = 1.0L / newest.inverse_curvature;
    const SpectralReal scale = yy > 1e-30L
        ? std::clamp(sy / yy, 1e-8L, 1e8L)
        : 1.0L;
    for (ComplexVector& value : direction) {
        for (SpectralComplex& component : value) {
            component *= scale;
        }
    }
    for (std::size_t index = 0; index < history.size(); ++index) {
        const LbfgsPair& pair = history[index];
        const SpectralReal beta = pair.inverse_curvature *
            increment_inner_product(pair.gradient_delta, direction);
        add_scaled_increment(
            direction, pair.state_delta, alpha[index] - beta);
    }
    return direction;
}

TriadPartition objective_partition(const std::string& objective) {
    if (objective == "critical-local-integral") {
        return TriadPartition::local;
    }
    if (objective == "critical-nonlocal-integral") {
        return TriadPartition::nonlocal;
    }
    return TriadPartition::all;
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
    const TriadPartition partition = objective_partition(options.objective);
    const SpectralReal initial_q =
        objective_.evaluate(state).energy_level_quantity;
    SpectralReal previous_integrand =
        objective_.evaluate(state, partition).critical_integrand;
    SpectralReal critical_integral = 0.0L;
    SpectralReal maximum_q = initial_q;
    for (int step = 0; step < options.trajectory_steps; ++step) {
        dynamics_.rk4_step(state, options.viscosity, options.time_step);
        const StaticObjective sample = objective_.evaluate(state, partition);
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
    if (options.objective == "critical-integral" ||
        options.objective == "critical-local-integral" ||
        options.objective == "critical-nonlocal-integral") {
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
    if (options.method != "steepest" && options.method != "lbfgs") {
        throw std::invalid_argument(
            "gradient-search method must be steepest or lbfgs");
    }
    if (options.lbfgs_history < 1 || options.lbfgs_history > 64) {
        throw std::invalid_argument(
            "L-BFGS history must be between 1 and 64");
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
    std::vector<LbfgsPair> lbfgs_history;
    SpectralState previous_state;
    SpectralIncrement previous_gradient;
    bool has_previous_point = false;

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
        } else if (options.objective == "critical-integral" ||
                   options.objective == "critical-local-integral" ||
                   options.objective == "critical-nonlocal-integral") {
            trajectory = adjoint_.critical_integral_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps,
                objective_partition(options.objective));
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
        SpectralIncrement projected_gradient = trajectory.initial_gradient;
        const SpectralReal gradient_norm = project_to_search_tangent(
            projected_gradient, result.state, sobolev);
        const SpectralReal sobolev_value = sobolev.value(result.state);
        result.final_projected_gradient_norm = gradient_norm;
        record.projected_gradient_norm = gradient_norm;
        record.sobolev_value = sobolev_value;
        if (!(gradient_norm > 1e-24L) || !std::isfinite(gradient_norm)) {
            record.objective_after = result.objective;
            result.trace.push_back(record);
            break;
        }

        if (options.method == "lbfgs" && has_previous_point) {
            SpectralIncrement state_delta = increment_difference(
                result.state.velocity, previous_state.velocity);
            project_to_search_tangent(state_delta, result.state, sobolev);
            SpectralIncrement transported_previous_gradient =
                previous_gradient;
            project_to_search_tangent(
                transported_previous_gradient, result.state, sobolev);
            SpectralIncrement gradient_delta = increment_difference(
                transported_previous_gradient, projected_gradient);
            project_to_search_tangent(
                gradient_delta, result.state, sobolev);
            const SpectralReal state_delta_norm = increment_norm(state_delta);
            const SpectralReal gradient_delta_norm =
                increment_norm(gradient_delta);
            const SpectralReal curvature = increment_inner_product(
                state_delta, gradient_delta);
            if (state_delta_norm > 0.0L && gradient_delta_norm > 0.0L &&
                curvature > 1e-12L * state_delta_norm *
                                gradient_delta_norm) {
                if (static_cast<int>(lbfgs_history.size()) ==
                    options.lbfgs_history) {
                    lbfgs_history.erase(lbfgs_history.begin());
                }
                lbfgs_history.push_back(LbfgsPair{
                    std::move(state_delta), std::move(gradient_delta),
                    1.0L / curvature});
            }
        }

        SpectralIncrement direction = options.method == "lbfgs"
            ? lbfgs_ascent_direction(projected_gradient, lbfgs_history)
            : projected_gradient;
        bool using_quasi_newton =
            options.method == "lbfgs" && !lbfgs_history.empty();
        SpectralReal direction_norm = project_to_search_tangent(
            direction, result.state, sobolev);
        const SpectralReal directional_derivative =
            increment_inner_product(projected_gradient, direction);
        if (!(direction_norm > 1e-24L) ||
            !(directional_derivative >
              1e-12L * gradient_norm * direction_norm)) {
            direction = projected_gradient;
            direction_norm = gradient_norm;
            lbfgs_history.clear();
            using_quasi_newton = false;
        }
        for (ComplexVector& value : direction) {
            for (SpectralComplex& component : value) {
                component /= direction_norm;
            }
        }

        auto line_search = [&](const SpectralIncrement& search_direction) {
            SpectralReal trial_step = next_step;
            for (int line_step = 0;
                 line_step < options.line_search_steps; ++line_step) {
                ++record.line_search_evaluations;
                SpectralState candidate = dynamics_.add_increment(
                    result.state, search_direction, trial_step);
                dynamics_.enforce_constraints(candidate);
                sobolev.retract(candidate, target_energy);
                const SpectralReal candidate_objective =
                    objective_value(candidate, options);
                ++result.trajectory_evaluations;
                const SpectralReal improvement_floor =
                    1e-13L * std::max(
                        1e-30L, std::abs(result.objective));
                if (std::isfinite(candidate_objective) &&
                    candidate_objective >
                        result.objective + improvement_floor) {
                    previous_state = result.state;
                    previous_gradient = projected_gradient;
                    has_previous_point = true;
                    result.state = std::move(candidate);
                    result.objective = candidate_objective;
                    ++result.accepted_steps;
                    next_step = std::min(
                        options.initial_step, 1.5L * trial_step);
                    record.accepted_step = trial_step;
                    record.accepted = true;
                    return true;
                }
                trial_step *= 0.5L;
            }
            return false;
        };

        bool accepted = line_search(direction);
        if (!accepted && using_quasi_newton) {
            direction = projected_gradient;
            for (ComplexVector& value : direction) {
                for (SpectralComplex& component : value) {
                    component /= gradient_norm;
                }
            }
            lbfgs_history.clear();
            record.used_steepest_fallback = true;
            accepted = line_search(direction);
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
