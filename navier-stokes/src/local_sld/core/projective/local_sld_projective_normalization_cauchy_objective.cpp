#include "local_sld_projective_normalization_cauchy_objective.hpp"

#include "local_quartic_closure_objective.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

struct ShapeRegion {
    std::vector<std::size_t> indices;
    std::vector<Shape> shapes;
};

struct CauchyGraph {
    ShapeRegion selected;
    ShapeRegion tail;
    SpectralIncrement au;
    SpectralIncrement aau;
    SpectralIncrement selected_b;
    SpectralIncrement selected_ab;
    SpectralIncrement tail_b;
    SpectralIncrement tail_ab;
    SpectralReal full_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_h1_norm2 = 0.0L;
    SpectralReal tail_h2_norm2 = 0.0L;
};

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective normalization Cauchy pairing layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective normalization Cauchy Laplacian layout mismatch");
    }
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
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "projective normalization Cauchy gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] +=
                factor * source[mode][coordinate];
        }
    }
}

CauchyGraph build_graph(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads) {
    CauchyGraph graph;
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    graph.selected.indices.reserve(groups.size());
    graph.selected.shapes.reserve(groups.size());
    for (std::size_t index = 0; index < groups.size(); ++index) {
        graph.selected.indices.push_back(index);
        graph.selected.shapes.push_back(
            groups[index].primitive_squared_lengths);
        if (groups[index].primitive_squared_lengths[2] >
            core_maximum_height) {
            graph.tail.indices.push_back(index);
            graph.tail.shapes.push_back(
                groups[index].primitive_squared_lengths);
        }
    }
    graph.au = laplacian_weight(state, state.velocity);
    graph.aau = laplacian_weight(state, graph.au);
    graph.selected_b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, graph.selected.indices,
            state.velocity, state.velocity, threads);
    graph.tail_b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, graph.tail.indices,
            state.velocity, state.velocity, threads);
    graph.selected_ab = laplacian_weight(state, graph.selected_b);
    graph.tail_ab = laplacian_weight(state, graph.tail_b);
    graph.full_stretching = LocalQuarticClosureObjective(
        dynamics, TriadPartition::local).evaluate(state).signed_stretching;
    graph.enstrophy = pairing(state.velocity, graph.au);
    graph.palinstrophy = pairing(graph.au, graph.au);
    graph.selected_h1_norm2 = pairing(
        graph.selected_b, graph.selected_ab);
    graph.tail_h2_norm2 = pairing(graph.tail_ab, graph.tail_ab);
    return graph;
}

}  // namespace

LocalSldProjectiveNormalizationCauchyObjective::
LocalSldProjectiveNormalizationCauchyObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads)
    : dynamics_(dynamics),
      selection_(selection),
      core_maximum_height_(core_maximum_height),
      threads_(threads) {
    if (core_maximum_height < 1) {
        throw std::invalid_argument(
            "projective normalization Cauchy height must be positive");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective normalization Cauchy threads must be 1..256");
    }
}

LocalSldProjectiveNormalizationCauchyObjectiveValue
LocalSldProjectiveNormalizationCauchyObjective::evaluate(
    const SpectralState& state) const {
    const CauchyGraph graph = build_graph(
        dynamics_, state, selection_, core_maximum_height_, threads_);
    LocalSldProjectiveNormalizationCauchyObjectiveValue result;
    result.core_maximum_height = core_maximum_height_;
    result.selected_shape_count = graph.selected.indices.size();
    result.tail_shape_count = graph.tail.indices.size();
    result.full_stretching = graph.full_stretching;
    result.enstrophy = graph.enstrophy;
    result.palinstrophy = graph.palinstrophy;
    result.selected_aggregate_h1_norm2 = graph.selected_h1_norm2;
    result.tail_aggregate_h2_norm2 = graph.tail_h2_norm2;
    const SpectralReal denominator =
        graph.enstrophy * graph.enstrophy * graph.enstrophy *
        graph.palinstrophy * graph.palinstrophy *
        graph.palinstrophy * graph.palinstrophy * graph.palinstrophy;
    if (denominator > 1e-60L) {
        result.squared_cauchy_bound_power_one =
            2.25L * graph.full_stretching * graph.full_stretching *
            graph.selected_h1_norm2 * graph.tail_h2_norm2 /
            denominator;
        result.cauchy_bound_power_one = std::sqrt(std::max(
            result.squared_cauchy_bound_power_one, 0.0L));
        result.finite = std::isfinite(
            result.squared_cauchy_bound_power_one);
    }
    return result;
}

SpectralIncrement
LocalSldProjectiveNormalizationCauchyObjective::gradient(
    const SpectralState& state) const {
    const CauchyGraph graph = build_graph(
        dynamics_, state, selection_, core_maximum_height_, threads_);
    SpectralIncrement result(state.waves.size());
    const SpectralReal denominator =
        graph.enstrophy * graph.enstrophy * graph.enstrophy *
        graph.palinstrophy * graph.palinstrophy *
        graph.palinstrophy * graph.palinstrophy * graph.palinstrophy;
    if (graph.selected.indices.empty() || graph.tail.indices.empty() ||
        !(denominator > 1e-60L) || graph.full_stretching == 0.0L ||
        !(graph.selected_h1_norm2 > 1e-60L) ||
        !(graph.tail_h2_norm2 > 1e-60L)) {
        return result;
    }
    const SpectralReal objective =
        2.25L * graph.full_stretching * graph.full_stretching *
        graph.selected_h1_norm2 * graph.tail_h2_norm2 / denominator;
    const auto& selected_aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, graph.selected.shapes);
    const auto& tail_aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, graph.tail.shapes);

    SpectralIncrement selected_h1_cotangent = graph.selected_ab;
    scale(selected_h1_cotangent, 2.0L);
    const SpectralIncrement selected_h1_gradient =
        ProjectiveAdvectionDecomposition::vjp(
            state, selected_aggregate.front(),
            selected_h1_cotangent, threads_);
    SpectralIncrement tail_h2_cotangent = laplacian_weight(
        state, graph.tail_ab);
    scale(tail_h2_cotangent, 2.0L);
    const SpectralIncrement tail_h2_gradient =
        ProjectiveAdvectionDecomposition::vjp(
            state, tail_aggregate.front(),
            tail_h2_cotangent, threads_);

    result = LocalQuarticClosureObjective(
        dynamics_, TriadPartition::local).signed_stretching_gradient(state);
    scale(result, 2.0L * objective / graph.full_stretching);
    add_scaled(
        result, selected_h1_gradient,
        objective / graph.selected_h1_norm2);
    add_scaled(
        result, tail_h2_gradient,
        objective / graph.tail_h2_norm2);
    add_scaled(result, graph.au, -6.0L * objective / graph.enstrophy);
    add_scaled(result, graph.aau, -10.0L * objective / graph.palinstrophy);
    return result;
}

}  // namespace lemma
