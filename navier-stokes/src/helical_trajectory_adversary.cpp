#include "helical_trajectory_adversary.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal inner_product(const SpectralIncrement& left,
                           const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal project_to_energy_sphere(
    SpectralIncrement& gradient, const SpectralState& state) {
    const SpectralReal energy = SpectralStateOps::energy(state);
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

}  // namespace

HelicalTrajectoryAdversaryResult HelicalTrajectoryAdversary::maximize(
    const SpectralState& initial,
    const HelicalTrajectoryAdversaryOptions& options,
    const HelicalSectorAdjoint& adjoint) {
    if (options.selection.sector_mask == 0U || options.iterations < 1 ||
        options.line_search_steps < 1 || options.trajectory_steps < 1 ||
        !(options.initial_step > 0.0L) || !(options.viscosity > 0.0L) ||
        !(options.time_step > 0.0L)) {
        throw std::invalid_argument("invalid helical trajectory options");
    }
    HelicalTrajectoryAdversaryResult result;
    result.state = initial;
    const SpectralReal energy = SpectralStateOps::energy(initial);
    result.objective = adjoint.critical_integral(
        result.state, options.viscosity, options.time_step,
        options.trajectory_steps, options.selection);
    result.initial_objective = result.objective;
    ++result.evaluations;
    for (int iteration = 0; iteration < options.iterations; ++iteration) {
        const HelicalSectorTrajectoryGradient gradient_result =
            adjoint.critical_integral_gradient(
                result.state, options.viscosity, options.time_step,
                options.trajectory_steps, options.selection);
        SpectralIncrement direction = gradient_result.initial_gradient;
        const SpectralReal gradient_norm = project_to_energy_sphere(
            direction, result.state);
        HelicalTrajectoryAdversaryTraceRow row;
        row.iteration = iteration;
        row.objective = result.objective;
        row.projected_gradient_norm = gradient_norm;
        if (!(gradient_norm > 0.0L)) {
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
            const SpectralReal candidate_objective = adjoint.critical_integral(
                candidate, options.viscosity, options.time_step,
                options.trajectory_steps, options.selection);
            ++result.evaluations;
            if (candidate_objective > result.objective) {
                result.state = std::move(candidate);
                result.objective = candidate_objective;
                row.accepted = true;
                row.accepted_step = step;
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

}  // namespace lemma
