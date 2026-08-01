#include "gradient_adversary.hpp"

#include <algorithm>
#include <cmath>
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

}  // namespace

GradientAdversary::GradientAdversary(const SpectralDynamics& dynamics,
                                     const SpectralObjective& objective,
                                     const SpectralAdjoint& adjoint)
    : dynamics_(dynamics), objective_(objective), adjoint_(adjoint) {}

SpectralReal GradientAdversary::maximum_q(
    const SpectralState& initial,
    const GradientSearchOptions& options) const {
    SpectralState state = initial;
    SpectralReal result =
        objective_.evaluate(state).energy_level_quantity;
    for (int step = 0; step < options.trajectory_steps; ++step) {
        dynamics_.rk4_step(state, options.viscosity, options.time_step);
        result = std::max(
            result, objective_.evaluate(state).energy_level_quantity);
    }
    return result;
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
    result.state = initial;
    const SpectralReal target_energy = SpectralStateOps::energy(initial);
    result.objective = maximum_q(result.state, options);
    result.initial_objective = result.objective;
    ++result.trajectory_evaluations;
    SpectralReal next_step = options.initial_step;

    for (int iteration = 0; iteration < options.iterations; ++iteration) {
        const QTrajectoryGradient trajectory = adjoint_.maximum_q_gradient(
            result.state, options.viscosity, options.time_step,
            options.trajectory_steps);
        ++result.trajectory_evaluations;
        ++result.iterations;
        result.objective = trajectory.objective_value;
        result.objective_step = trajectory.objective_step;
        SpectralIncrement direction = trajectory.initial_gradient;
        const SpectralReal gradient_norm =
            project_to_energy_sphere(direction, result.state);
        result.final_projected_gradient_norm = gradient_norm;
        if (!(gradient_norm > 1e-24L) || !std::isfinite(gradient_norm)) {
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
            SpectralState candidate = dynamics_.add_increment(
                result.state, direction, trial_step);
            SpectralStateOps::normalize_energy(candidate, target_energy);
            const SpectralReal candidate_objective = maximum_q(candidate, options);
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
                accepted = true;
                break;
            }
            trial_step *= 0.5L;
        }
        if (!accepted) {
            break;
        }
    }
    return result;
}

}  // namespace lemma
