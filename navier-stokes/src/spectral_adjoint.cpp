#include "spectral_adjoint.hpp"

#include <stdexcept>
#include <cmath>
#include <utility>
#include <vector>

namespace lemma {
namespace {

std::vector<SpectralState> build_checkpoints(
    const SpectralDynamics& dynamics, const SpectralState& initial,
    SpectralReal viscosity, SpectralReal time_step, int steps) {
    if (!(viscosity > 0.0L)) {
        throw std::invalid_argument("adjoint viscosity must be positive");
    }
    if (!(time_step > 0.0L)) {
        throw std::invalid_argument("adjoint time step must be positive");
    }
    if (steps < 0) {
        throw std::invalid_argument("adjoint step count cannot be negative");
    }
    std::vector<SpectralState> checkpoints;
    checkpoints.reserve(static_cast<std::size_t>(steps) + 1);
    checkpoints.push_back(initial);
    for (int step = 0; step < steps; ++step) {
        SpectralState next = checkpoints.back();
        dynamics.rk4_step(next, viscosity, time_step);
        checkpoints.push_back(std::move(next));
    }
    return checkpoints;
}

}  // namespace

SpectralAdjoint::SpectralAdjoint(const SpectralDynamics& dynamics,
                                 const SpectralObjective& objective)
    : dynamics_(dynamics), objective_(objective) {}

QTrajectoryGradient SpectralAdjoint::terminal_q_gradient(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps) const {
    const std::vector<SpectralState> checkpoints = build_checkpoints(
        dynamics_, initial, viscosity, time_step, steps);
    return reverse_from_step(checkpoints, viscosity, time_step, steps);
}

QTrajectoryGradient SpectralAdjoint::maximum_q_gradient(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps) const {
    const std::vector<SpectralState> checkpoints = build_checkpoints(
        dynamics_, initial, viscosity, time_step, steps);
    int maximum_step = 0;
    SpectralReal maximum_q =
        objective_.evaluate(checkpoints.front()).energy_level_quantity;
    for (int step = 1; step <= steps; ++step) {
        const SpectralReal q = objective_
            .evaluate(checkpoints[static_cast<std::size_t>(step)])
            .energy_level_quantity;
        if (q > maximum_q) {
            maximum_q = q;
            maximum_step = step;
        }
    }
    return reverse_from_step(checkpoints, viscosity, time_step, maximum_step);
}

QTrajectoryGradient SpectralAdjoint::q_gain_gradient(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps) const {
    const std::vector<SpectralState> checkpoints = build_checkpoints(
        dynamics_, initial, viscosity, time_step, steps);
    QTrajectoryGradient result = reverse_from_step(
        checkpoints, viscosity, time_step, steps);
    const SpectralReal initial_q =
        objective_.evaluate(initial).energy_level_quantity;
    const SpectralReal terminal_q = result.objective_value;
    if (!(initial_q > 1e-30L) || !(terminal_q > 1e-30L)) {
        throw std::runtime_error(
            "q-gain gradient requires positive endpoint Q values");
    }
    const SpectralIncrement initial_q_gradient =
        objective_.energy_level_gradient(initial);
    for (std::size_t mode = 0; mode < result.initial_gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            result.initial_gradient[mode][component] =
                result.initial_gradient[mode][component] / terminal_q -
                initial_q_gradient[mode][component] / initial_q;
        }
    }
    result.objective_value = std::log(terminal_q / initial_q);
    return result;
}

QTrajectoryGradient SpectralAdjoint::q_increase_gradient(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps) const {
    const std::vector<SpectralState> checkpoints = build_checkpoints(
        dynamics_, initial, viscosity, time_step, steps);
    QTrajectoryGradient result = reverse_from_step(
        checkpoints, viscosity, time_step, steps);
    const SpectralReal initial_q =
        objective_.evaluate(initial).energy_level_quantity;
    const SpectralReal terminal_q = result.objective_value;
    const SpectralIncrement initial_q_gradient =
        objective_.energy_level_gradient(initial);
    for (std::size_t mode = 0; mode < result.initial_gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            result.initial_gradient[mode][component] -=
                initial_q_gradient[mode][component];
        }
    }
    result.objective_value = terminal_q - initial_q;
    return result;
}

QTrajectoryGradient SpectralAdjoint::critical_integral_gradient(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps) const {
    const std::vector<SpectralState> checkpoints = build_checkpoints(
        dynamics_, initial, viscosity, time_step, steps);
    QTrajectoryGradient result;
    result.objective_step = steps;
    result.total_steps = steps;
    result.checkpoint_count = checkpoints.size();
    result.initial_gradient = SpectralIncrement(initial.waves.size());
    if (steps == 0) {
        return result;
    }
    auto add_scaled = [](SpectralIncrement& target,
                         const SpectralIncrement& source,
                         SpectralReal scale) {
        for (std::size_t mode = 0; mode < target.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                target[mode][component] +=
                    scale * source[mode][component];
            }
        }
    };
    for (int step = 0; step <= steps; ++step) {
        const SpectralReal weight =
            (step == 0 || step == steps) ? 0.5L * time_step : time_step;
        result.objective_value +=
            weight * objective_
                         .evaluate(checkpoints[static_cast<std::size_t>(step)])
                         .critical_integrand;
    }
    result.initial_gradient = objective_.critical_integrand_gradient(
        checkpoints.back());
    for (ComplexVector& value : result.initial_gradient) {
        for (SpectralComplex& component : value) {
            component *= 0.5L * time_step;
        }
    }
    for (int step = steps - 1; step >= 0; --step) {
        result.initial_gradient = dynamics_.rk4_vjp(
            checkpoints[static_cast<std::size_t>(step)],
            result.initial_gradient, viscosity, time_step);
        const SpectralReal weight =
            step == 0 ? 0.5L * time_step : time_step;
        const SpectralIncrement source =
            objective_.critical_integrand_gradient(
                checkpoints[static_cast<std::size_t>(step)]);
        add_scaled(result.initial_gradient, source, weight);
    }
    return result;
}

QTrajectoryGradient SpectralAdjoint::reverse_from_step(
    const std::vector<SpectralState>& checkpoints,
    SpectralReal viscosity, SpectralReal time_step,
    int objective_step) const {
    if (objective_step < 0 ||
        static_cast<std::size_t>(objective_step) >= checkpoints.size()) {
        throw std::out_of_range("adjoint objective checkpoint is out of range");
    }
    QTrajectoryGradient result;
    result.objective_step = objective_step;
    result.total_steps = static_cast<int>(checkpoints.size()) - 1;
    result.checkpoint_count = checkpoints.size();
    const SpectralState& objective_state =
        checkpoints[static_cast<std::size_t>(objective_step)];
    result.objective_value =
        objective_.evaluate(objective_state).energy_level_quantity;
    result.initial_gradient =
        objective_.energy_level_gradient(objective_state);
    for (int step = objective_step - 1; step >= 0; --step) {
        result.initial_gradient = dynamics_.rk4_vjp(
            checkpoints[static_cast<std::size_t>(step)],
            result.initial_gradient, viscosity, time_step);
    }
    return result;
}

}  // namespace lemma
