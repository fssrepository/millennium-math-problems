#include "local_sld_projective_height_stretching_objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

struct ShellGroups {
    std::vector<std::size_t> indices;
    std::vector<Shape> shapes;
};

ShellGroups shell_groups(
    const std::vector<ProjectiveInteractionGroup>& groups,
    SpectralInteger lower,
    SpectralInteger upper) {
    ShellGroups result;
    for (std::size_t index = 0; index < groups.size(); ++index) {
        const SpectralInteger height =
            groups[index].primitive_squared_lengths[2];
        if (height > lower && height <= upper) {
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
            "projective height-stretching pairing layout mismatch");
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

LocalSldProjectiveHeightStretchingObjective::
LocalSldProjectiveHeightStretchingObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection,
    SpectralInteger core_maximum_height,
    int threads)
    : dynamics_(dynamics),
      selection_(selection),
      core_maximum_height_(core_maximum_height),
      shell_maximum_height_(2 * core_maximum_height),
      threads_(threads) {
    if (core_maximum_height < 1) {
        throw std::invalid_argument(
            "projective height-stretching core height must be positive");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height-stretching threads must be 1..256");
    }
}

LocalSldProjectiveHeightStretchingObjectiveValue
LocalSldProjectiveHeightStretchingObjective::evaluate(
    const SpectralState& state) const {
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    const ShellGroups shell = shell_groups(
        groups, core_maximum_height_, shell_maximum_height_);
    LocalSldProjectiveHeightStretchingObjectiveValue result;
    result.core_maximum_height = core_maximum_height_;
    result.shell_maximum_height = shell_maximum_height_;
    result.shell_shape_count = shell.indices.size();
    if (shell.indices.empty()) {
        return result;
    }
    const SpectralIncrement b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, shell.indices,
            state.velocity, state.velocity, threads_);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const SpectralIncrement ab = laplacian_weight(state, b);
    result.signed_shell_stretching = pairing(au, b);
    result.enstrophy = pairing(state.velocity, au);
    result.aggregate_h1_norm2 = pairing(b, ab);
    result.square_function_h1_norm2 =
        ProjectiveAdvectionDecomposition::square_function_norms(
            state, groups, shell.indices, threads_).h1_norm2;
    if (!(result.enstrophy > 1e-60L) ||
        !(result.square_function_h1_norm2 > 1e-60L)) {
        return result;
    }
    result.stretching_aware_h1_ratio =
        result.signed_shell_stretching *
        result.signed_shell_stretching /
        (result.enstrophy * result.square_function_h1_norm2);
    if (result.aggregate_h1_norm2 > 1e-60L) {
        result.h1_synthesis_ratio =
            result.aggregate_h1_norm2 /
            result.square_function_h1_norm2;
        result.stretching_h1_alignment_squared =
            result.signed_shell_stretching *
            result.signed_shell_stretching /
            (result.enstrophy * result.aggregate_h1_norm2);
        const SpectralReal reconstructed =
            result.h1_synthesis_ratio *
            result.stretching_h1_alignment_squared;
        result.product_reconstruction_error = std::abs(
            reconstructed - result.stretching_aware_h1_ratio) /
            std::max({std::abs(reconstructed),
                      std::abs(result.stretching_aware_h1_ratio),
                      1e-30L});
    }
    result.finite = std::isfinite(result.stretching_aware_h1_ratio) &&
        std::isfinite(result.product_reconstruction_error);
    return result;
}

SpectralIncrement
LocalSldProjectiveHeightStretchingObjective::gradient(
    const SpectralState& state) const {
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection_);
    const ShellGroups shell = shell_groups(
        groups, core_maximum_height_, shell_maximum_height_);
    SpectralIncrement result(state.waves.size());
    if (shell.indices.empty()) {
        return result;
    }
    const auto& aggregate =
        ProjectiveAdvectionDecomposition::aggregate_family(
            state, selection_, shell.shapes);
    const SpectralIncrement b =
        ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, shell.indices,
            state.velocity, state.velocity, threads_);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const SpectralReal stretching = pairing(au, b);
    const SpectralReal enstrophy = pairing(state.velocity, au);
    ProjectiveSquareFunctionMoment square_function =
        ProjectiveAdvectionDecomposition::h1_square_function(
            state, groups, shell.indices, true);
    if (!(enstrophy > 1e-60L) ||
        !(square_function.norm2 > 1e-60L)) {
        return result;
    }
    SpectralIncrement stretching_gradient =
        laplacian_weight(state, b);
    add_scaled(
        stretching_gradient,
        ProjectiveAdvectionDecomposition::vjp(
            state, aggregate.front(), au),
        1.0L);
    const SpectralReal denominator =
        enstrophy * square_function.norm2;
    const SpectralReal objective =
        stretching * stretching / denominator;
    result = stretching_gradient;
    scale(result, 2.0L * stretching / denominator);
    add_scaled(result, au, -2.0L * objective / enstrophy);
    add_scaled(
        result, square_function.gradient,
        -objective / square_function.norm2);
    return result;
}

}  // namespace lemma
