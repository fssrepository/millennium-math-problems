#include "local_sld_projective_height_dynamic_ratio_objective.hpp"

#include "projective_height_envelope_kernel.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

void scale(SpectralIncrement& value, SpectralReal factor) {
    for (ComplexVector& mode : value) {
        for (SpectralComplex& component : mode) {
            component *= factor;
        }
    }
}

void add_scaled(
    SpectralIncrement& target,
    const SpectralIncrement& source,
    SpectralReal factor) {
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] +=
                factor * source[mode][coordinate];
        }
    }
}

}  // namespace

LocalSldProjectiveHeightDynamicRatioObjective::
LocalSldProjectiveHeightDynamicRatioObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    int threads)
    : dynamics_(dynamics), selection_(selection), threads_(threads) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height dynamic-ratio threads must be 1..256");
    }
}

LocalSldProjectiveHeightDynamicRatioObjectiveValue
LocalSldProjectiveHeightDynamicRatioObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    LocalSldProjectiveHeightDynamicRatioObjectiveValue result;
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L)) {
        return result;
    }
    const ProjectiveHeightEnvelopeMoment envelope =
        ProjectiveHeightEnvelopeKernel::evaluate(
            state, selection_, selected.enstrophy,
            selected.palinstrophy, false, threads_, true, true);
    result.height_shell_count = envelope.height_shell_count;
    result.dynamic_paired_bracket_envelope =
        envelope.absolute_component_envelope;
    result.outer_h1_sum = envelope.outer_h1_sum;
    if (!(result.outer_h1_sum > 1e-30L)) {
        return result;
    }
    result.coercivity_ratio =
        result.dynamic_paired_bracket_envelope / result.outer_h1_sum;
    result.squared_coercivity_ratio =
        result.coercivity_ratio * result.coercivity_ratio;
    result.finite = std::isfinite(result.squared_coercivity_ratio);
    return result;
}

SpectralIncrement
LocalSldProjectiveHeightDynamicRatioObjective::gradient(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    SpectralIncrement result(state.waves.size());
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L)) {
        return result;
    }
    const ProjectiveHeightEnvelopeMoment envelope =
        ProjectiveHeightEnvelopeKernel::evaluate(
            state, selection_, selected.enstrophy,
            selected.palinstrophy, true, threads_, true, true);
    if (!((envelope.outer_h1_sum > 1e-30L))) {
        return result;
    }
    const SpectralReal ratio =
        envelope.absolute_component_envelope / envelope.outer_h1_sum;
    result = envelope.gradient;
    scale(result, 1.0L / envelope.outer_h1_sum);
    add_scaled(
        result, envelope.outer_gradient,
        -envelope.absolute_component_envelope /
            (envelope.outer_h1_sum * envelope.outer_h1_sum));
    scale(result, 2.0L * ratio);
    return result;
}

}  // namespace lemma
