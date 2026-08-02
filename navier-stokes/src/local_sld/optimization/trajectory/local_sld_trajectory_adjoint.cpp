#include "local_sld_trajectory_adjoint.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lemma {
namespace {

void validate(SpectralReal viscosity,
              SpectralReal time_step,
              int steps) {
    if (!(viscosity > 0.0L) || !(time_step > 0.0L) || steps < 0) {
        throw std::invalid_argument(
            "frozen SLD trajectory requires positive viscosity/time step and nonnegative steps");
    }
}

void add_scaled(SpectralIncrement& target,
                const SpectralIncrement& source,
                SpectralReal factor) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "frozen SLD trajectory gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += factor * source[mode][component];
        }
    }
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source,
    int power) {
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        SpectralReal weight = 1.0L;
        for (int exponent = 0; exponent < power; ++exponent) {
            weight *= wave2;
        }
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

LocalSldTrajectoryValue evaluate_terminal(
    const LocalQuarticClosureObjective& objective,
    const SpectralState& initial,
    const SpectralState& terminal,
    int steps) {
    const LocalQuarticClosureObjectiveValue initial_value =
        objective.evaluate(initial);
    const LocalQuarticClosureObjectiveValue terminal_value =
        objective.evaluate(terminal);
    LocalSldTrajectoryValue result;
    result.initial_frequency = std::sqrt(
        initial_value.enstrophy / initial_value.energy);
    result.initial_ep_shift =
        initial_value.energy * initial_value.palinstrophy;
    result.terminal_ratio = objective.frozen_signed_local_sld_ratio(
        terminal, result.initial_frequency, result.initial_ep_shift);
    const SpectralReal s2 = terminal_value.signed_stretching *
        terminal_value.signed_stretching;
    const SpectralReal s4 = s2 * s2;
    const SpectralReal z = terminal_value.enstrophy;
    const SpectralReal p = terminal_value.palinstrophy;
    const SpectralReal p2 = p * p;
    const SpectralReal p4 = p2 * p2;
    const SpectralReal first = s4 * z * z * p;
    const SpectralReal shifted = result.initial_ep_shift *
        z * z * z * p4;
    if (first + shifted > 0.0L) {
        result.terminal_shift_fraction = shifted / (first + shifted);
    }
    result.steps = steps;
    result.finite = initial_value.finite && terminal_value.finite &&
        std::isfinite(result.terminal_ratio) &&
        result.initial_frequency > 0.0L &&
        result.initial_ep_shift > 0.0L;
    return result;
}

}  // namespace

LocalSldTrajectoryAdjoint::LocalSldTrajectoryAdjoint(
    const SpectralDynamics& dynamics,
    TriadSelection selection)
    : dynamics_(dynamics), selection_(selection) {}

LocalSldTrajectoryValue LocalSldTrajectoryAdjoint::terminal_value(
    const SpectralState& initial,
    SpectralReal viscosity,
    SpectralReal time_step,
    int steps) const {
    validate(viscosity, time_step, steps);
    const LocalQuarticClosureObjective objective(dynamics_, selection_);
    SpectralState terminal = initial;
    for (int step = 0; step < steps; ++step) {
        dynamics_.rk4_step(terminal, viscosity, time_step);
    }
    return evaluate_terminal(objective, initial, terminal, steps);
}

QTrajectoryGradient LocalSldTrajectoryAdjoint::terminal_gradient(
    const SpectralState& initial,
    SpectralReal viscosity,
    SpectralReal time_step,
    int steps) const {
    validate(viscosity, time_step, steps);
    const LocalQuarticClosureObjective objective(dynamics_, selection_);
    std::vector<SpectralState> checkpoints;
    checkpoints.reserve(static_cast<std::size_t>(steps) + 1);
    checkpoints.push_back(initial);
    for (int step = 0; step < steps; ++step) {
        SpectralState next = checkpoints.back();
        dynamics_.rk4_step(next, viscosity, time_step);
        checkpoints.push_back(std::move(next));
    }
    const LocalSldTrajectoryValue value = evaluate_terminal(
        objective, checkpoints.front(), checkpoints.back(), steps);
    QTrajectoryGradient result;
    result.objective_value = value.terminal_ratio;
    result.objective_step = steps;
    result.total_steps = steps;
    result.checkpoint_count = checkpoints.size();
    result.initial_gradient =
        objective.frozen_signed_local_sld_ratio_gradient(
            checkpoints.back(), value.initial_frequency,
            value.initial_ep_shift);
    for (int step = steps - 1; step >= 0; --step) {
        result.initial_gradient = dynamics_.rk4_vjp(
            checkpoints[static_cast<std::size_t>(step)],
            result.initial_gradient, viscosity, time_step);
    }

    const LocalQuarticClosureObjectiveValue initial_value =
        objective.evaluate(initial);
    const SpectralIncrement au = laplacian_weight(
        initial, initial.velocity, 1);
    const SpectralIncrement a2u = laplacian_weight(
        initial, initial.velocity, 2);
    add_scaled(
        result.initial_gradient, au,
        -value.terminal_ratio / initial_value.enstrophy);
    add_scaled(
        result.initial_gradient, initial.velocity,
        value.terminal_ratio / initial_value.energy);
    const SpectralReal shift_weight =
        -2.0L * value.terminal_ratio *
        value.terminal_shift_fraction;
    add_scaled(
        result.initial_gradient, initial.velocity,
        shift_weight / initial_value.energy);
    add_scaled(
        result.initial_gradient, a2u,
        shift_weight / initial_value.palinstrophy);
    return result;
}

LocalSldTrajectoryValue LocalSldTrajectoryAdjoint::maximum_value(
    const SpectralState& initial,
    SpectralReal viscosity,
    SpectralReal time_step,
    int steps) const {
    validate(viscosity, time_step, steps);
    const LocalQuarticClosureObjective objective(dynamics_, selection_);
    SpectralState state = initial;
    LocalSldTrajectoryValue best = evaluate_terminal(
        objective, initial, state, 0);
    for (int step = 1; step <= steps; ++step) {
        dynamics_.rk4_step(state, viscosity, time_step);
        const LocalSldTrajectoryValue candidate = evaluate_terminal(
            objective, initial, state, step);
        if (candidate.terminal_ratio > best.terminal_ratio) {
            best = candidate;
        }
    }
    return best;
}

QTrajectoryGradient LocalSldTrajectoryAdjoint::maximum_gradient(
    const SpectralState& initial,
    SpectralReal viscosity,
    SpectralReal time_step,
    int steps) const {
    validate(viscosity, time_step, steps);
    const LocalQuarticClosureObjective objective(dynamics_, selection_);
    std::vector<SpectralState> checkpoints;
    checkpoints.reserve(static_cast<std::size_t>(steps) + 1);
    checkpoints.push_back(initial);
    LocalSldTrajectoryValue best = evaluate_terminal(
        objective, initial, initial, 0);
    int best_step = 0;
    for (int step = 1; step <= steps; ++step) {
        SpectralState next = checkpoints.back();
        dynamics_.rk4_step(next, viscosity, time_step);
        checkpoints.push_back(std::move(next));
        const LocalSldTrajectoryValue candidate = evaluate_terminal(
            objective, initial, checkpoints.back(), step);
        if (candidate.terminal_ratio > best.terminal_ratio) {
            best = candidate;
            best_step = step;
        }
    }
    QTrajectoryGradient result;
    result.objective_value = best.terminal_ratio;
    result.objective_step = best_step;
    result.total_steps = steps;
    result.checkpoint_count = checkpoints.size();
    result.initial_gradient =
        objective.frozen_signed_local_sld_ratio_gradient(
            checkpoints[static_cast<std::size_t>(best_step)],
            best.initial_frequency, best.initial_ep_shift);
    for (int step = best_step - 1; step >= 0; --step) {
        result.initial_gradient = dynamics_.rk4_vjp(
            checkpoints[static_cast<std::size_t>(step)],
            result.initial_gradient, viscosity, time_step);
    }

    const LocalQuarticClosureObjectiveValue initial_value =
        objective.evaluate(initial);
    const SpectralIncrement au = laplacian_weight(
        initial, initial.velocity, 1);
    const SpectralIncrement a2u = laplacian_weight(
        initial, initial.velocity, 2);
    add_scaled(
        result.initial_gradient, au,
        -best.terminal_ratio / initial_value.enstrophy);
    add_scaled(
        result.initial_gradient, initial.velocity,
        best.terminal_ratio / initial_value.energy);
    const SpectralReal shift_weight =
        -2.0L * best.terminal_ratio * best.terminal_shift_fraction;
    add_scaled(
        result.initial_gradient, initial.velocity,
        shift_weight / initial_value.energy);
    add_scaled(
        result.initial_gradient, a2u,
        shift_weight / initial_value.palinstrophy);
    return result;
}

}  // namespace lemma
