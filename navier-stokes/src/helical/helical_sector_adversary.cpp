#include "helical_sector_adversary.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal inner_product(const SpectralIncrement& left,
                           const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument("helical adversary layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal project_to_energy_sphere(
    SpectralIncrement& gradient, const SpectralState& state) {
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

}  // namespace

HelicalSectorAdversaryResult HelicalSectorAdversary::maximize(
    const SpectralState& initial,
    const HelicalSectorAdversaryOptions& options) {
    if (options.selection.sector_mask == 0U || options.iterations < 1 ||
        options.line_search_steps < 1 || !(options.initial_step > 0.0L)) {
        throw std::invalid_argument("invalid helical adversary options");
    }
    HelicalSectorAdversaryResult result;
    result.state = initial;
    const SpectralReal energy = SpectralStateOps::energy(initial);
    result.objective = HelicalSectorObjective::evaluate(
        result.state, options.selection).critical_integrand;
    result.initial_objective = result.objective;
    ++result.evaluations;
    for (int iteration = 0; iteration < options.iterations; ++iteration) {
        SpectralIncrement direction =
            HelicalSectorObjective::critical_integrand_gradient(
                result.state, options.selection);
        const SpectralReal gradient_norm = project_to_energy_sphere(
            direction, result.state);
        HelicalSectorAdversaryTraceRow row;
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
            const SpectralReal candidate_objective =
                HelicalSectorObjective::evaluate(
                    candidate, options.selection).critical_integrand;
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

}  // namespace lemma
