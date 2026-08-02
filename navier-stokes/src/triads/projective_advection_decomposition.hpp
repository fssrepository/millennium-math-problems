#pragma once

#include "spectral_state.hpp"
#include "triad_partition.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace lemma {

struct ProjectiveInteractionGroup {
    std::array<SpectralInteger, 3> primitive_squared_lengths{};
    std::vector<InteractionIndex> interactions;
};

struct ProjectiveSquareFunctionMoment {
    SpectralReal norm2 = 0.0L;
    SpectralIncrement gradient;
};

struct ProjectiveSquareFunctionNorms {
    SpectralReal l2_norm2 = 0.0L;
    SpectralReal h1_norm2 = 0.0L;
    SpectralReal h2_norm2 = 0.0L;
};

struct ProjectiveBilinearCotangents {
    SpectralIncrement advecting;
    SpectralIncrement advected;
};

class ProjectiveAdvectionDecomposition {
public:
    [[nodiscard]] static const std::vector<ProjectiveInteractionGroup>& group(
        const SpectralState& state,
        TriadSelection selection);

    [[nodiscard]] static const std::vector<ProjectiveInteractionGroup>&
    aggregate_family(
        const SpectralState& state,
        TriadSelection selection,
        std::vector<std::array<SpectralInteger, 3>>
            primitive_squared_lengths);

    [[nodiscard]] static SpectralIncrement evaluate(
        const SpectralState& state,
        const ProjectiveInteractionGroup& group);

    [[nodiscard]] static SpectralIncrement evaluate_bilinear(
        const SpectralState& state,
        const ProjectiveInteractionGroup& group,
        const SpectralIncrement& advecting,
        const SpectralIncrement& advected,
        int threads = 1);

    [[nodiscard]] static SpectralIncrement evaluate_bilinear_sum(
        const SpectralState& state,
        const std::vector<ProjectiveInteractionGroup>& groups,
        const std::vector<std::size_t>& group_indices,
        const SpectralIncrement& advecting,
        const SpectralIncrement& advected,
        int threads = 1);

    [[nodiscard]] static SpectralIncrement vjp(
        const SpectralState& state,
        const ProjectiveInteractionGroup& group,
        const SpectralIncrement& output_cotangent,
        int threads = 1);

    [[nodiscard]] static SpectralIncrement vjp_sum(
        const SpectralState& state,
        const std::vector<ProjectiveInteractionGroup>& groups,
        const std::vector<std::size_t>& group_indices,
        const SpectralIncrement& output_cotangent,
        int threads = 1);

    [[nodiscard]] static ProjectiveBilinearCotangents bilinear_vjp(
        const SpectralState& state,
        const ProjectiveInteractionGroup& group,
        const SpectralIncrement& advecting,
        const SpectralIncrement& advected,
        const SpectralIncrement& output_cotangent,
        int threads = 1);

    [[nodiscard]] static ProjectiveSquareFunctionMoment square_function(
        const SpectralState& state,
        const std::vector<ProjectiveInteractionGroup>& groups,
        bool compute_gradient);

    [[nodiscard]] static ProjectiveSquareFunctionNorms
    square_function_norms(
        const SpectralState& state,
        const std::vector<ProjectiveInteractionGroup>& groups,
        const std::vector<std::size_t>& group_indices,
        int threads = 1);

    [[nodiscard]] static ProjectiveSquareFunctionMoment
    h1_square_function(
        const SpectralState& state,
        const std::vector<ProjectiveInteractionGroup>& groups,
        const std::vector<std::size_t>& group_indices,
        bool compute_gradient);
};

}  // namespace lemma
