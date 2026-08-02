#include "local_sld_projective_stretching_objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective stretching pairing layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal norm2(const SpectralIncrement& value) {
    return pairing(value, value);
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective stretching Laplacian layout mismatch");
    }
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= wave2;
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
            "projective stretching gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += factor * source[mode][component];
        }
    }
}

struct ObjectiveGraph {
    std::size_t projective_shape_count = 0;
    SpectralIncrement au;
    SpectralIncrement coherent;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal coherent_norm2 = 0.0L;
    SpectralReal square_function_norm2 = 0.0L;
};

ObjectiveGraph build_graph(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    TriadSelection selection) {
    ObjectiveGraph graph;
    const std::vector<ProjectiveInteractionGroup>& groups =
        ProjectiveAdvectionDecomposition::group(
        state, selection);
    graph.projective_shape_count = groups.size();
    graph.au = laplacian_weight(state, state.velocity);
    graph.coherent = dynamics.advection_direct_partition(
        state, selection);
    graph.stretching = pairing(graph.au, graph.coherent);
    graph.palinstrophy = norm2(graph.au);
    graph.coherent_norm2 = norm2(graph.coherent);
    graph.square_function_norm2 =
        ProjectiveAdvectionDecomposition::square_function(
            state, groups, false).norm2;
    return graph;
}

}  // namespace

LocalSldProjectiveStretchingObjective::
LocalSldProjectiveStretchingObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection)
    : dynamics_(dynamics), selection_(selection) {}

LocalSldProjectiveStretchingObjectiveValue
LocalSldProjectiveStretchingObjective::evaluate(
    const SpectralState& state) const {
    const ObjectiveGraph graph = build_graph(
        dynamics_, state, selection_);
    LocalSldProjectiveStretchingObjectiveValue result;
    result.projective_shape_count = graph.projective_shape_count;
    result.signed_selected_stretching = graph.stretching;
    result.palinstrophy = graph.palinstrophy;
    result.coherent_norm2 = graph.coherent_norm2;
    result.square_function_norm2 = graph.square_function_norm2;
    if (!(graph.palinstrophy > 1e-60L) ||
        !(graph.square_function_norm2 > 1e-60L)) {
        return result;
    }
    result.stretching_aware_synthesis_ratio =
        graph.stretching * graph.stretching /
        (graph.palinstrophy * graph.square_function_norm2);
    if (graph.coherent_norm2 > 1e-60L) {
        result.coherent_synthesis_ratio =
            graph.coherent_norm2 / graph.square_function_norm2;
        result.stretching_alignment_squared =
            graph.stretching * graph.stretching /
            (graph.palinstrophy * graph.coherent_norm2);
        const SpectralReal reconstructed =
            result.coherent_synthesis_ratio *
            result.stretching_alignment_squared;
        result.product_reconstruction_error = std::abs(
            reconstructed - result.stretching_aware_synthesis_ratio) /
            std::max({std::abs(reconstructed),
                      std::abs(result.stretching_aware_synthesis_ratio),
                      1e-30L});
    }
    result.finite = std::isfinite(
        result.stretching_aware_synthesis_ratio) &&
        std::isfinite(result.product_reconstruction_error);
    return result;
}

SpectralIncrement LocalSldProjectiveStretchingObjective::gradient(
    const SpectralState& state) const {
    const ObjectiveGraph graph = build_graph(
        dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    if (!(graph.palinstrophy > 1e-60L) ||
        !(graph.square_function_norm2 > 1e-60L)) {
        return result;
    }

    SpectralIncrement stretching_gradient = laplacian_weight(
        state, graph.coherent);
    add_scaled(
        stretching_gradient,
        dynamics_.advection_vjp_direct_partition(
            state, graph.au, selection_),
        1.0L);

    const std::vector<ProjectiveInteractionGroup>& groups =
        ProjectiveAdvectionDecomposition::group(state, selection_);
    ProjectiveSquareFunctionMoment square_function =
        ProjectiveAdvectionDecomposition::square_function(
            state, groups, true);

    const SpectralReal denominator =
        graph.palinstrophy * graph.square_function_norm2;
    const SpectralReal objective =
        graph.stretching * graph.stretching / denominator;
    result = stretching_gradient;
    scale(result, 2.0L * graph.stretching / denominator);
    add_scaled(
        result,
        laplacian_weight(state, graph.au),
        -2.0L * objective / graph.palinstrophy);
    add_scaled(
        result,
        square_function.gradient,
        -objective / graph.square_function_norm2);
    return result;
}

}  // namespace lemma
