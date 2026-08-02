#include "local_sld_projective_normalization_schur_objective.hpp"

#include "local_quartic_closure_objective.hpp"
#include "projective_advection_decomposition.hpp"
#include "projective_height_shell_partition.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lemma {
namespace {

struct SchurShellGraph {
    int shell = 0;
    std::vector<std::size_t> group_indices;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralReal h2_norm2 = 0.0L;
};

struct SchurGraph {
    std::vector<std::size_t> selected_indices;
    SpectralIncrement au;
    SpectralIncrement aau;
    SpectralIncrement selected_b;
    SpectralIncrement selected_ab;
    std::vector<SchurShellGraph> shells;
    std::size_t row_index = 0;
    SpectralReal full_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_h1_norm2 = 0.0L;
    SpectralReal diagonal_tail_h2_norm2 = 0.0L;
    SpectralReal row_sum = 0.0L;
};

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "projective normalization Schur pairing layout mismatch");
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
            "projective normalization Schur Laplacian layout mismatch");
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
            "projective normalization Schur gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += factor * source[mode][component];
        }
    }
}

SchurGraph build_graph(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int row_shell,
    int threads) {
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    SchurGraph graph;
    graph.selected_indices.reserve(groups.size());
    std::map<int, SchurShellGraph> tail;
    for (std::size_t index = 0; index < groups.size(); ++index) {
        graph.selected_indices.push_back(index);
        const SpectralInteger height =
            groups[index].primitive_squared_lengths[2];
        if (height <= core_maximum_height) {
            continue;
        }
        const int shell = ProjectiveHeightShellPartition::shell_index(
            height);
        auto& target = tail[shell];
        target.shell = shell;
        target.group_indices.push_back(index);
    }
    if (tail.find(row_shell) == tail.end()) {
        throw std::invalid_argument(
            "projective normalization Schur row is absent from the tail");
    }
    graph.au = laplacian_weight(state, state.velocity);
    graph.aau = laplacian_weight(state, graph.au);
    graph.selected_b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, graph.selected_indices,
            state.velocity, state.velocity, threads);
    graph.selected_ab = laplacian_weight(state, graph.selected_b);
    graph.selected_h1_norm2 = pairing(
        graph.selected_b, graph.selected_ab);
    graph.shells.reserve(tail.size());
    for (auto& [shell_index, shell] : tail) {
        static_cast<void>(shell_index);
        shell.b = ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, shell.group_indices,
            state.velocity, state.velocity, threads);
        shell.ab = laplacian_weight(state, shell.b);
        shell.h2_norm2 = pairing(shell.ab, shell.ab);
        graph.diagonal_tail_h2_norm2 += shell.h2_norm2;
        if (shell.shell == row_shell) {
            graph.row_index = graph.shells.size();
        }
        graph.shells.push_back(std::move(shell));
    }
    const auto& row = graph.shells[graph.row_index];
    if (row.shell != row_shell || !(row.h2_norm2 > 1e-60L)) {
        throw std::invalid_argument(
            "projective normalization Schur row has zero output");
    }
    graph.row_sum = 1.0L;
    for (std::size_t index = 0; index < graph.shells.size(); ++index) {
        if (index == graph.row_index ||
            !(graph.shells[index].h2_norm2 > 1e-60L)) {
            continue;
        }
        const SpectralReal denominator = std::sqrt(
            row.h2_norm2 * graph.shells[index].h2_norm2);
        graph.row_sum += std::abs(pairing(
            row.ab, graph.shells[index].ab)) / denominator;
    }
    graph.full_stretching = LocalQuarticClosureObjective(
        dynamics, TriadPartition::local).evaluate(state).signed_stretching;
    graph.enstrophy = pairing(state.velocity, graph.au);
    graph.palinstrophy = pairing(graph.au, graph.au);
    return graph;
}

SpectralReal common_factor(const SchurGraph& graph) {
    const SpectralReal denominator =
        graph.enstrophy * graph.enstrophy * graph.enstrophy *
        graph.palinstrophy * graph.palinstrophy * graph.palinstrophy *
        graph.palinstrophy * graph.palinstrophy;
    if (!(denominator > 1e-60L)) {
        return 0.0L;
    }
    return 2.25L * graph.full_stretching * graph.full_stretching *
        graph.selected_h1_norm2 / denominator;
}

}  // namespace

LocalSldProjectiveNormalizationSchurObjective::
LocalSldProjectiveNormalizationSchurObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int row_shell,
    int threads)
    : dynamics_(dynamics), selection_(selection),
      core_maximum_height_(core_maximum_height), row_shell_(row_shell),
      threads_(threads) {
    if (core_maximum_height < 1 || row_shell < 0) {
        throw std::invalid_argument(
            "projective normalization Schur requires positive height and nonnegative row");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective normalization Schur threads must be 1..256");
    }
}

LocalSldProjectiveNormalizationSchurObjectiveValue
LocalSldProjectiveNormalizationSchurObjective::evaluate(
    const SpectralState& state) const {
    const SchurGraph graph = build_graph(
        dynamics_, state, selection_, core_maximum_height_,
        row_shell_, threads_);
    LocalSldProjectiveNormalizationSchurObjectiveValue result;
    result.core_maximum_height = core_maximum_height_;
    result.row_shell = row_shell_;
    result.selected_shape_count = graph.selected_indices.size();
    result.tail_shell_count = graph.shells.size();
    result.full_stretching = graph.full_stretching;
    result.enstrophy = graph.enstrophy;
    result.palinstrophy = graph.palinstrophy;
    result.selected_aggregate_h1_norm2 = graph.selected_h1_norm2;
    result.diagonal_tail_h2_norm2 = graph.diagonal_tail_h2_norm2;
    result.normalized_absolute_gram_row_sum = graph.row_sum;
    result.normalization_common_factor = common_factor(graph);
    result.schur_squared_majorant = result.normalization_common_factor *
        graph.row_sum * graph.diagonal_tail_h2_norm2;
    result.height_half_compensated_schur_squared_majorant = std::sqrt(
        static_cast<SpectralReal>(core_maximum_height_)) *
        result.schur_squared_majorant;
    result.finite =
        result.normalization_common_factor >= 0.0L &&
        std::isfinite(
            result.height_half_compensated_schur_squared_majorant);
    return result;
}

SpectralIncrement
LocalSldProjectiveNormalizationSchurObjective::gradient(
    const SpectralState& state) const {
    const SchurGraph graph = build_graph(
        dynamics_, state, selection_, core_maximum_height_,
        row_shell_, threads_);
    SpectralIncrement result(state.waves.size());
    const SpectralReal common = common_factor(graph);
    const SpectralReal objective = std::sqrt(
        static_cast<SpectralReal>(core_maximum_height_)) * common *
        graph.row_sum * graph.diagonal_tail_h2_norm2;
    if (!(objective > 0.0L) || graph.full_stretching == 0.0L ||
        !(graph.selected_h1_norm2 > 1e-60L) ||
        !(graph.diagonal_tail_h2_norm2 > 1e-60L) ||
        !(graph.row_sum > 0.0L)) {
        return result;
    }

    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    SpectralIncrement selected_h1_cotangent = graph.selected_ab;
    scale(selected_h1_cotangent, 2.0L);
    const SpectralIncrement selected_h1_gradient =
        ProjectiveAdvectionDecomposition::vjp_sum(
            state, groups, graph.selected_indices,
            selected_h1_cotangent, threads_);
    result = LocalQuarticClosureObjective(
        dynamics_, TriadPartition::local).signed_stretching_gradient(state);
    scale(result, 2.0L * objective / graph.full_stretching);
    add_scaled(
        result, selected_h1_gradient,
        objective / graph.selected_h1_norm2);
    add_scaled(result, graph.au, -6.0L * objective / graph.enstrophy);
    add_scaled(result, graph.aau, -10.0L * objective / graph.palinstrophy);

    std::vector<SpectralIncrement> row_cotangents;
    row_cotangents.reserve(graph.shells.size());
    for (std::size_t index = 0; index < graph.shells.size(); ++index) {
        row_cotangents.emplace_back(state.waves.size());
    }
    const auto& row = graph.shells[graph.row_index];
    for (std::size_t index = 0; index < graph.shells.size(); ++index) {
        if (index == graph.row_index ||
            !(graph.shells[index].h2_norm2 > 1e-60L)) {
            continue;
        }
        const auto& other = graph.shells[index];
        const SpectralReal gram = pairing(row.ab, other.ab);
        if (gram == 0.0L) {
            continue;
        }
        const SpectralReal sign = gram > 0.0L ? 1.0L : -1.0L;
        const SpectralReal denominator = std::sqrt(
            row.h2_norm2 * other.h2_norm2);
        const SpectralReal correlation = std::abs(gram) / denominator;
        add_scaled(
            row_cotangents[graph.row_index], other.ab,
            sign / denominator);
        add_scaled(
            row_cotangents[graph.row_index], row.ab,
            -correlation / row.h2_norm2);
        add_scaled(
            row_cotangents[index], row.ab,
            sign / denominator);
        add_scaled(
            row_cotangents[index], other.ab,
            -correlation / other.h2_norm2);
    }

    for (std::size_t index = 0; index < graph.shells.size(); ++index) {
        const auto& shell = graph.shells[index];
        SpectralIncrement output_cotangent = laplacian_weight(
            state, shell.ab);
        scale(
            output_cotangent,
            2.0L * objective / graph.diagonal_tail_h2_norm2);
        SpectralIncrement row_output_cotangent = laplacian_weight(
            state, row_cotangents[index]);
        add_scaled(
            output_cotangent, row_output_cotangent,
            objective / graph.row_sum);
        const SpectralIncrement shell_gradient =
            ProjectiveAdvectionDecomposition::vjp_sum(
                state, groups, shell.group_indices,
                output_cotangent, threads_);
        add_scaled(result, shell_gradient, 1.0L);
    }
    return result;
}

}  // namespace lemma
