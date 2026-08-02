#include "projective_advection_decomposition.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

Shape primitive_shape(
    const SpectralState& state,
    InteractionIndex interaction) {
    const auto [p, q, target] = interaction;
    Shape signature{
        norm_squared(state.waves[p]),
        norm_squared(state.waves[q]),
        norm_squared(state.waves[target])};
    std::sort(signature.begin(), signature.end());
    const SpectralInteger divisor = std::gcd(
        signature[0], std::gcd(signature[1], signature[2]));
    if (divisor <= 0) {
        throw std::invalid_argument(
            "projective advection requires nonzero wave vectors");
    }
    return {signature[0] / divisor,
            signature[1] / divisor,
            signature[2] / divisor};
}

void require_layout(
    const SpectralState& state,
    const SpectralIncrement& increment) {
    if (increment.size() != state.waves.size()) {
        throw std::invalid_argument(
            "projective advection layout mismatch");
    }
}

void project(
    SpectralIncrement& increment,
    const SpectralState& state) {
    for (std::size_t mode = 0; mode < increment.size(); ++mode) {
        increment[mode] = project_divergence_free(
            state.waves[mode], increment[mode]);
    }
}

}  // namespace

std::vector<ProjectiveInteractionGroup>
ProjectiveAdvectionDecomposition::group(
    const SpectralState& state,
    TriadSelection selection) {
    std::map<Shape, std::vector<InteractionIndex>> grouped;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (TriadPartitioner::includes(state, interaction, selection)) {
            grouped[primitive_shape(state, interaction)].push_back(
                interaction);
        }
    }
    std::vector<ProjectiveInteractionGroup> result;
    result.reserve(grouped.size());
    for (auto& [shape, interactions] : grouped) {
        result.push_back({shape, std::move(interactions)});
    }
    return result;
}

SpectralIncrement ProjectiveAdvectionDecomposition::evaluate(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group) {
    SpectralIncrement result(state.waves.size());
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const SpectralComplex coefficient = imaginary_unit *
            wave_dot(state.waves[q], state.velocity[p]);
        for (std::size_t component = 0; component < 3; ++component) {
            result[target][component] +=
                coefficient * state.velocity[q][component];
        }
    }
    project(result, state);
    return result;
}

SpectralIncrement ProjectiveAdvectionDecomposition::vjp(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SpectralIncrement& output_cotangent) {
    require_layout(state, output_cotangent);
    SpectralIncrement cotangent = output_cotangent;
    project(cotangent, state);
    SpectralIncrement result(state.waves.size());
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const ComplexVector& target_cotangent = cotangent[target];
        const SpectralComplex first_coefficient =
            minus_imaginary_unit *
            dot_hermitian(state.velocity[q], target_cotangent);
        for (std::size_t component = 0; component < 3; ++component) {
            const SpectralReal wave_component =
                static_cast<SpectralReal>(
                    component == 0   ? state.waves[q].x
                    : component == 1 ? state.waves[q].y
                                     : state.waves[q].z);
            result[p][component] +=
                wave_component * first_coefficient;
        }
        const SpectralComplex second_coefficient =
            minus_imaginary_unit * std::conj(
                wave_dot(state.waves[q], state.velocity[p]));
        for (std::size_t component = 0; component < 3; ++component) {
            result[q][component] +=
                second_coefficient * target_cotangent[component];
        }
    }
    project(result, state);
    return result;
}

}  // namespace lemma
