#include "local_sld_projective_coherence_objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

SpectralReal norm2(const SpectralIncrement& value) {
    SpectralReal result = 0.0L;
    for (const ComplexVector& mode : value) {
        for (const SpectralComplex component : mode) {
            result += std::norm(component);
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
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "projective coherence gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += factor * source[mode][component];
        }
    }
}

}  // namespace

LocalSldProjectiveCoherenceObjective::
LocalSldProjectiveCoherenceObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection)
    : dynamics_(dynamics), selection_(selection) {}

LocalSldProjectiveCoherenceObjectiveValue
LocalSldProjectiveCoherenceObjective::evaluate(
    const SpectralState& state) const {
    LocalSldProjectiveCoherenceObjectiveValue result;
    const std::vector<ProjectiveInteractionGroup>& groups =
        ProjectiveAdvectionDecomposition::group(state, selection_);
    result.projective_shape_count = groups.size();
    const SpectralIncrement coherent =
        dynamics_.advection_direct_partition(state, selection_);
    result.coherent_norm2 = norm2(coherent);
    result.square_function_norm2 =
        ProjectiveAdvectionDecomposition::square_function(
            state, groups, false).norm2;
    if (result.square_function_norm2 > 1e-60L) {
        result.synthesis_ratio = result.coherent_norm2 /
            result.square_function_norm2;
        result.synthesis_amplification = std::sqrt(std::max(
            0.0L, result.synthesis_ratio));
        result.finite = std::isfinite(result.synthesis_ratio);
    }
    return result;
}

SpectralIncrement LocalSldProjectiveCoherenceObjective::gradient(
    const SpectralState& state) const {
    const std::vector<ProjectiveInteractionGroup>& groups =
        ProjectiveAdvectionDecomposition::group(state, selection_);
    const SpectralIncrement coherent =
        dynamics_.advection_direct_partition(state, selection_);
    const SpectralReal numerator = norm2(coherent);
    ProjectiveSquareFunctionMoment square_function =
        ProjectiveAdvectionDecomposition::square_function(
            state, groups, true);
    const SpectralReal denominator = square_function.norm2;
    if (!(denominator > 1e-60L)) {
        return SpectralIncrement(state.waves.size());
    }
    SpectralIncrement numerator_cotangent = coherent;
    scale(numerator_cotangent, 2.0L);
    SpectralIncrement result =
        dynamics_.advection_vjp_direct_partition(
            state, numerator_cotangent, selection_);
    scale(result, 1.0L / denominator);
    add_scaled(
        result, square_function.gradient,
        -numerator / (denominator * denominator));
    return result;
}

}  // namespace lemma
