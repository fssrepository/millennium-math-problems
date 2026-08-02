#include "local_sld_projective_core_tail_alignment_objective.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

struct RegionGroups {
    std::vector<std::size_t> indices;
    std::vector<Shape> shapes;
};

RegionGroups region_groups(
    const std::vector<ProjectiveInteractionGroup>& groups,
    SpectralInteger core_maximum_height,
    LocalSldProjectiveHeightRegion region) {
    RegionGroups result;
    for (std::size_t index = 0; index < groups.size(); ++index) {
        const SpectralInteger height =
            groups[index].primitive_squared_lengths[2];
        const bool in_core = height <= core_maximum_height;
        if ((region == LocalSldProjectiveHeightRegion::core && in_core) ||
            (region == LocalSldProjectiveHeightRegion::tail && !in_core)) {
            result.indices.push_back(index);
            result.shapes.push_back(
                groups[index].primitive_squared_lengths);
        }
    }
    return result;
}

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective core-tail alignment pairing layout mismatch");
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
            "projective core-tail alignment Laplacian layout mismatch");
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
            "projective core-tail alignment gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] +=
                factor * source[mode][coordinate];
        }
    }
}

struct ObjectiveGraph {
    RegionGroups region;
    SpectralIncrement au;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralReal stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal aggregate_h1_norm2 = 0.0L;
};

ObjectiveGraph build_graph(
    const SpectralState& state,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    LocalSldProjectiveHeightRegion region,
    int threads) {
    ObjectiveGraph graph;
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    graph.region = region_groups(
        groups, core_maximum_height, region);
    graph.au = laplacian_weight(state, state.velocity);
    graph.b = ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
        state, groups, graph.region.indices,
        state.velocity, state.velocity, threads);
    graph.ab = laplacian_weight(state, graph.b);
    graph.stretching = pairing(graph.au, graph.b);
    graph.enstrophy = pairing(state.velocity, graph.au);
    graph.aggregate_h1_norm2 = pairing(graph.b, graph.ab);
    return graph;
}

}  // namespace

LocalSldProjectiveCoreTailAlignmentObjective::
LocalSldProjectiveCoreTailAlignmentObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    LocalSldProjectiveHeightRegion region,
    int threads)
    : dynamics_(dynamics),
      selection_(selection),
      core_maximum_height_(core_maximum_height),
      region_(region),
      threads_(threads) {
    if (core_maximum_height < 1) {
        throw std::invalid_argument(
            "projective core-tail alignment height must be positive");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective core-tail alignment threads must be 1..256");
    }
}

LocalSldProjectiveCoreTailAlignmentObjectiveValue
LocalSldProjectiveCoreTailAlignmentObjective::evaluate(
    const SpectralState& state) const {
    const ObjectiveGraph graph = build_graph(
        state, selection_, core_maximum_height_, region_, threads_);
    LocalSldProjectiveCoreTailAlignmentObjectiveValue result;
    result.core_maximum_height = core_maximum_height_;
    result.region = region_;
    result.projective_shape_count = graph.region.indices.size();
    result.signed_stretching = graph.stretching;
    result.enstrophy = graph.enstrophy;
    result.aggregate_h1_norm2 = graph.aggregate_h1_norm2;
    const SpectralReal denominator =
        graph.enstrophy * graph.aggregate_h1_norm2;
    if (denominator > 1e-60L) {
        result.stretching_h1_alignment_squared =
            graph.stretching * graph.stretching / denominator;
        result.finite = std::isfinite(
            result.stretching_h1_alignment_squared);
    }
    return result;
}

SpectralIncrement
LocalSldProjectiveCoreTailAlignmentObjective::gradient(
    const SpectralState& state) const {
    const ObjectiveGraph graph = build_graph(
        state, selection_, core_maximum_height_, region_, threads_);
    SpectralIncrement result(state.waves.size());
    const SpectralReal denominator =
        graph.enstrophy * graph.aggregate_h1_norm2;
    if (graph.region.indices.empty() || !(denominator > 1e-60L)) {
        return result;
    }
    const auto& aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, graph.region.shapes);

    SpectralIncrement stretching_gradient =
        laplacian_weight(state, graph.b);
    add_scaled(
        stretching_gradient,
        ProjectiveAdvectionDecomposition::vjp(
            state, aggregate.front(), graph.au, threads_),
        1.0L);

    SpectralIncrement h1_cotangent = graph.ab;
    scale(h1_cotangent, 2.0L);
    const SpectralIncrement h1_gradient =
        ProjectiveAdvectionDecomposition::vjp(
            state, aggregate.front(), h1_cotangent, threads_);

    const SpectralReal objective =
        graph.stretching * graph.stretching / denominator;
    result = stretching_gradient;
    scale(
        result,
        2.0L * graph.stretching / denominator);
    add_scaled(
        result, graph.au,
        -2.0L * objective / graph.enstrophy);
    add_scaled(
        result, h1_gradient,
        -objective / graph.aggregate_h1_norm2);
    return result;
}

}  // namespace lemma
