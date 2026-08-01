#include "helical_sector_adjoint.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace lemma {
namespace {

void add_scaled(SpectralIncrement& target, const SpectralIncrement& source,
                SpectralReal scale) {
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += scale * source[mode][component];
        }
    }
}

}  // namespace

SpectralReal HelicalSectorAdjoint::critical_integral(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps,
    HelicalSectorSelection selection) const {
    if (!(viscosity > 0.0L) || !(time_step > 0.0L) || steps < 1 ||
        selection.sector_mask == 0U) {
        throw std::invalid_argument("invalid helical integral parameters");
    }
    SpectralState state = initial;
    SpectralReal integral = 0.5L * time_step *
        HelicalSectorObjective::evaluate(state, selection).critical_integrand;
    for (int step = 0; step < steps; ++step) {
        dynamics_.rk4_step(state, viscosity, time_step);
        const SpectralReal weight = step + 1 == steps
            ? 0.5L * time_step
            : time_step;
        integral += weight *
            HelicalSectorObjective::evaluate(state, selection)
                .critical_integrand;
    }
    return integral;
}

HelicalSectorTrajectoryGradient
HelicalSectorAdjoint::critical_integral_gradient(
    const SpectralState& initial, SpectralReal viscosity,
    SpectralReal time_step, int steps,
    HelicalSectorSelection selection) const {
    if (!(viscosity > 0.0L) || !(time_step > 0.0L) || steps < 1 ||
        selection.sector_mask == 0U) {
        throw std::invalid_argument("invalid helical adjoint parameters");
    }
    std::vector<SpectralState> checkpoints;
    checkpoints.reserve(static_cast<std::size_t>(steps) + 1);
    checkpoints.push_back(initial);
    for (int step = 0; step < steps; ++step) {
        SpectralState next = checkpoints.back();
        dynamics_.rk4_step(next, viscosity, time_step);
        checkpoints.push_back(std::move(next));
    }

    HelicalSectorTrajectoryGradient result;
    result.total_steps = steps;
    result.checkpoint_count = checkpoints.size();
    for (int step = 0; step <= steps; ++step) {
        const SpectralReal weight =
            (step == 0 || step == steps) ? 0.5L * time_step : time_step;
        result.objective_value += weight *
            HelicalSectorObjective::evaluate(
                checkpoints[static_cast<std::size_t>(step)], selection)
                .critical_integrand;
    }

    result.initial_gradient =
        HelicalSectorObjective::critical_integrand_gradient(
            checkpoints.back(), selection);
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
        add_scaled(
            result.initial_gradient,
            HelicalSectorObjective::critical_integrand_gradient(
                checkpoints[static_cast<std::size_t>(step)], selection),
            weight);
    }
    return result;
}

}  // namespace lemma
