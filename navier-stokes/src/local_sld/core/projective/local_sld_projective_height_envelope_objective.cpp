#include "local_sld_projective_height_envelope_objective.hpp"

#include "projective_height_envelope_kernel.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

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

SpectralReal sign(SpectralReal value) {
    return value > 0.0L ? 1.0L : (value < 0.0L ? -1.0L : 0.0L);
}

}  // namespace

LocalSldProjectiveHeightEnvelopeObjective::
LocalSldProjectiveHeightEnvelopeObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    int threads,
    bool pair_outer_and_advected_commutator)
    : dynamics_(dynamics),
      selection_(selection),
      threads_(threads),
      pair_outer_and_advected_commutator_(
          pair_outer_and_advected_commutator) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height envelope threads must be 1..256");
    }
}

LocalSldProjectiveHeightEnvelopeObjectiveValue
LocalSldProjectiveHeightEnvelopeObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics_, TriadPartition::local).evaluate(state);
    LocalSldProjectiveHeightEnvelopeObjectiveValue result;
    result.enstrophy = selected.enstrophy;
    result.palinstrophy = selected.palinstrophy;
    result.full_stretching = full.signed_stretching;
    result.pairs_outer_and_advected_commutator =
        pair_outer_and_advected_commutator_;
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    const ProjectiveHeightEnvelopeMoment envelope =
        ProjectiveHeightEnvelopeKernel::evaluate(
            state, selection_, result.enstrophy,
            result.palinstrophy, false, threads_,
            pair_outer_and_advected_commutator_);
    result.height_shell_count = envelope.height_shell_count;
    result.height_pair_count = envelope.height_pair_count;
    result.absolute_component_brackets =
        envelope.absolute_component_sums;
    result.absolute_component_bracket_envelope =
        envelope.absolute_component_envelope;
    const SpectralReal denominator =
        result.enstrophy * result.enstrophy *
        result.palinstrophy * result.palinstrophy;
    result.absolute_component_power_one_envelope =
        std::abs(result.full_stretching) *
        result.absolute_component_bracket_envelope / denominator;
    result.squared_component_power_one_envelope =
        result.absolute_component_power_one_envelope *
        result.absolute_component_power_one_envelope;
    result.finite = std::isfinite(
        result.squared_component_power_one_envelope);
    return result;
}

SpectralIncrement LocalSldProjectiveHeightEnvelopeObjective::gradient(
    const SpectralState& state) const {
    const LocalQuarticClosureObjective selected_objective(
        dynamics_, selection_);
    const LocalQuarticClosureObjective full_objective(
        dynamics_, TriadPartition::local);
    const LocalQuarticClosureObjectiveValue selected =
        selected_objective.evaluate(state);
    const LocalQuarticClosureObjectiveValue full =
        full_objective.evaluate(state);
    SpectralIncrement result(state.waves.size());
    if (!(selected.enstrophy > 0.0L) ||
        !(selected.palinstrophy > 0.0L)) {
        return result;
    }
    const ProjectiveHeightEnvelopeMoment envelope =
        ProjectiveHeightEnvelopeKernel::evaluate(
            state, selection_, selected.enstrophy,
            selected.palinstrophy, true, threads_,
            pair_outer_and_advected_commutator_);
    const SpectralReal denominator =
        selected.enstrophy * selected.enstrophy *
        selected.palinstrophy * selected.palinstrophy;
    const SpectralReal absolute_stretching =
        std::abs(full.signed_stretching);
    const SpectralReal power = absolute_stretching *
        envelope.absolute_component_envelope / denominator;
    result = envelope.gradient;
    scale(result, absolute_stretching / denominator);
    add_scaled(
        result,
        full_objective.signed_stretching_gradient(state),
        sign(full.signed_stretching) *
            envelope.absolute_component_envelope / denominator);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    add_scaled(result, au, -4.0L * power / selected.enstrophy);
    add_scaled(
        result, laplacian_weight(state, au),
        -4.0L * power / selected.palinstrophy);
    scale(result, 2.0L * power);
    return result;
}

}  // namespace lemma
