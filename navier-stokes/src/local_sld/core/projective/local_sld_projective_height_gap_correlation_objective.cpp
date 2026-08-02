#include "local_sld_projective_height_gap_correlation_objective.hpp"

#include "projective_height_shell_partition.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

struct GapGraph {
    const std::vector<ProjectiveInteractionGroup>* groups = nullptr;
    std::vector<std::size_t> first_group_indices;
    std::vector<std::size_t> second_group_indices;
    SpectralIncrement first_ab;
    SpectralIncrement second_ab;
    SpectralReal gram = 0.0L;
    SpectralReal first_norm2 = 0.0L;
    SpectralReal second_norm2 = 0.0L;
};

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "height-gap correlation pairing layout mismatch");
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

void add(
    SpectralIncrement& target,
    const SpectralIncrement& source) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "height-gap correlation gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            target[mode][coordinate] += source[mode][coordinate];
        }
    }
}

GapGraph build_graph(
    const SpectralState& state,
    TriadSelection selection,
    int first_shell,
    int second_shell,
    int threads) {
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    const auto partition = ProjectiveHeightShellPartition::build(groups);
    if (first_shell < 0 || second_shell < 0 ||
        static_cast<std::size_t>(first_shell) >= partition.size() ||
        static_cast<std::size_t>(second_shell) >= partition.size() ||
        partition[static_cast<std::size_t>(first_shell)]
            .group_indices.empty() ||
        partition[static_cast<std::size_t>(second_shell)]
            .group_indices.empty()) {
        throw std::invalid_argument(
            "height-gap correlation shell is absent from the cutoff layout");
    }
    const auto& first_partition =
        partition[static_cast<std::size_t>(first_shell)];
    const auto& second_partition =
        partition[static_cast<std::size_t>(second_shell)];
    const SpectralIncrement first_b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, first_partition.group_indices,
            state.velocity, state.velocity, threads);
    const SpectralIncrement second_b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, second_partition.group_indices,
            state.velocity, state.velocity, threads);
    GapGraph graph;
    graph.groups = &groups;
    graph.first_group_indices = first_partition.group_indices;
    graph.second_group_indices = second_partition.group_indices;
    graph.first_ab = laplacian_weight(state, first_b);
    graph.second_ab = laplacian_weight(state, second_b);
    graph.gram = pairing(graph.first_ab, graph.second_ab);
    graph.first_norm2 = pairing(graph.first_ab, graph.first_ab);
    graph.second_norm2 = pairing(graph.second_ab, graph.second_ab);
    return graph;
}

}  // namespace

LocalSldProjectiveHeightGapCorrelationObjective::
LocalSldProjectiveHeightGapCorrelationObjective(
    TriadSelection selection,
    int first_shell,
    int second_shell,
    int threads)
    : selection_(selection),
      first_shell_(first_shell),
      second_shell_(second_shell),
      threads_(threads) {
    if (first_shell < 0 || second_shell <= first_shell ||
        second_shell > 62) {
        throw std::invalid_argument(
            "height-gap correlation requires 0 <= first < second <= 62");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "height-gap correlation threads must be 1..256");
    }
}

LocalSldProjectiveHeightGapCorrelationObjectiveValue
LocalSldProjectiveHeightGapCorrelationObjective::evaluate(
    const SpectralState& state) const {
    const GapGraph graph = build_graph(
        state, selection_, first_shell_, second_shell_, threads_);
    LocalSldProjectiveHeightGapCorrelationObjectiveValue result;
    result.first_shell = first_shell_;
    result.second_shell = second_shell_;
    result.shell_gap = second_shell_ - first_shell_;
    result.gram_pairing = graph.gram;
    result.first_h2_norm2 = graph.first_norm2;
    result.second_h2_norm2 = graph.second_norm2;
    const SpectralReal denominator =
        graph.first_norm2 * graph.second_norm2;
    if (denominator > 1e-60L) {
        result.correlation_squared =
            graph.gram * graph.gram / denominator;
        result.half_decay_weighted_correlation = std::ldexp(
            std::sqrt(std::max(result.correlation_squared, 0.0L)),
            result.shell_gap);
        result.weighted_correlation_squared = std::ldexp(
            result.correlation_squared, 2 * result.shell_gap);
        result.finite = std::isfinite(
            result.weighted_correlation_squared);
    }
    return result;
}

SpectralIncrement
LocalSldProjectiveHeightGapCorrelationObjective::gradient(
    const SpectralState& state) const {
    const GapGraph graph = build_graph(
        state, selection_, first_shell_, second_shell_, threads_);
    SpectralIncrement result(state.waves.size());
    const SpectralReal denominator =
        graph.first_norm2 * graph.second_norm2;
    if (!(denominator > 1e-60L)) {
        return result;
    }
    const SpectralReal correlation_squared =
        graph.gram * graph.gram / denominator;
    const SpectralReal weight = std::ldexp(
        1.0L, 2 * (second_shell_ - first_shell_));
    SpectralIncrement first_ab_cotangent = graph.second_ab;
    scale(
        first_ab_cotangent,
        weight * 2.0L * graph.gram / denominator);
    SpectralIncrement first_self = graph.first_ab;
    scale(
        first_self,
        -weight * 2.0L * correlation_squared / graph.first_norm2);
    add(first_ab_cotangent, first_self);

    SpectralIncrement second_ab_cotangent = graph.first_ab;
    scale(
        second_ab_cotangent,
        weight * 2.0L * graph.gram / denominator);
    SpectralIncrement second_self = graph.second_ab;
    scale(
        second_self,
        -weight * 2.0L * correlation_squared / graph.second_norm2);
    add(second_ab_cotangent, second_self);

    const SpectralIncrement first_b_cotangent = laplacian_weight(
        state, first_ab_cotangent);
    const SpectralIncrement second_b_cotangent = laplacian_weight(
        state, second_ab_cotangent);
    result = ProjectiveAdvectionDecomposition::vjp_sum(
        state, *graph.groups, graph.first_group_indices,
        first_b_cotangent, threads_);
    add(
        result,
        ProjectiveAdvectionDecomposition::vjp_sum(
            state, *graph.groups, graph.second_group_indices,
            second_b_cotangent, threads_));
    return result;
}

}  // namespace lemma
