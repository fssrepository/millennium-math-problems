#include "projective_advection_decomposition.hpp"

#include <algorithm>
#include <compare>
#include <map>
#include <mutex>
#include <numeric>
#include <stdexcept>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

struct GroupCacheKey {
    std::vector<WaveVector> waves;
    int minimum_gap = 0;
    int maximum_gap = 0;
    TriadSelection::SignatureMode signature_mode =
        TriadSelection::SignatureMode::none;
    Shape signature{};
    SpectralInteger minimum_low_squared = 0;
    SpectralInteger maximum_low_squared_exclusive = 0;

    auto operator<=>(const GroupCacheKey&) const = default;
};

struct FamilyGroupCacheKey {
    GroupCacheKey group;
    std::vector<Shape> primitive_squared_lengths;

    auto operator<=>(const FamilyGroupCacheKey&) const = default;
};

GroupCacheKey cache_key(
    const SpectralState& state,
    TriadSelection selection) {
    return {
        state.waves,
        selection.minimum_gap(),
        selection.maximum_gap(),
        selection.signature_mode(),
        selection.squared_length_signature(),
        selection.minimum_low_squared(),
        selection.maximum_low_squared_exclusive()};
}

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

const std::vector<ProjectiveInteractionGroup>&
ProjectiveAdvectionDecomposition::group(
    const SpectralState& state,
    TriadSelection selection) {
    static std::map<
        GroupCacheKey, std::vector<ProjectiveInteractionGroup>> cache;
    static std::mutex cache_mutex;
    const GroupCacheKey key = cache_key(state, selection);
    const std::lock_guard<std::mutex> lock(cache_mutex);
    const auto existing = cache.find(key);
    if (existing != cache.end()) {
        return existing->second;
    }
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
    return cache.emplace(key, std::move(result)).first->second;
}

const std::vector<ProjectiveInteractionGroup>&
ProjectiveAdvectionDecomposition::aggregate_family(
    const SpectralState& state,
    TriadSelection selection,
    std::vector<Shape> primitive_squared_lengths) {
    for (Shape& shape : primitive_squared_lengths) {
        std::sort(shape.begin(), shape.end());
    }
    std::sort(
        primitive_squared_lengths.begin(),
        primitive_squared_lengths.end());
    if (std::adjacent_find(
            primitive_squared_lengths.begin(),
            primitive_squared_lengths.end()) !=
        primitive_squared_lengths.end()) {
        throw std::invalid_argument(
            "projective aggregate family contains duplicate shapes");
    }
    const auto& groups = group(state, selection);
    static std::map<
        FamilyGroupCacheKey,
        std::vector<ProjectiveInteractionGroup>> cache;
    static std::mutex cache_mutex;
    const FamilyGroupCacheKey key{
        cache_key(state, selection), primitive_squared_lengths};
    const std::lock_guard<std::mutex> lock(cache_mutex);
    const auto existing = cache.find(key);
    if (existing != cache.end()) {
        return existing->second;
    }
    std::vector<ProjectiveInteractionGroup> result(1);
    ProjectiveInteractionGroup& aggregate = result.front();
    std::size_t interaction_count = 0;
    for (const ProjectiveInteractionGroup& source : groups) {
        if (std::binary_search(
                primitive_squared_lengths.begin(),
                primitive_squared_lengths.end(),
                source.primitive_squared_lengths)) {
            interaction_count += source.interactions.size();
        }
    }
    aggregate.interactions.reserve(interaction_count);
    for (const ProjectiveInteractionGroup& source : groups) {
        if (std::binary_search(
                primitive_squared_lengths.begin(),
                primitive_squared_lengths.end(),
                source.primitive_squared_lengths)) {
            aggregate.interactions.insert(
                aggregate.interactions.end(),
                source.interactions.begin(),
                source.interactions.end());
        }
    }
    return cache.emplace(key, std::move(result)).first->second;
}

SpectralIncrement ProjectiveAdvectionDecomposition::evaluate(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group) {
    return evaluate_bilinear(
        state, group, state.velocity, state.velocity);
}

SpectralIncrement ProjectiveAdvectionDecomposition::evaluate_bilinear(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SpectralIncrement& advecting,
    const SpectralIncrement& advected) {
    require_layout(state, advecting);
    require_layout(state, advected);
    SpectralIncrement result(state.waves.size());
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const SpectralComplex coefficient = imaginary_unit *
            wave_dot(state.waves[q], advecting[p]);
        for (std::size_t component = 0; component < 3; ++component) {
            result[target][component] +=
                coefficient * advected[q][component];
        }
    }
    project(result, state);
    return result;
}

SpectralIncrement
ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    const std::vector<std::size_t>& group_indices,
    const SpectralIncrement& advecting,
    const SpectralIncrement& advected,
    int threads) {
    require_layout(state, advecting);
    require_layout(state, advected);
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective bilinear-sum threads must be 1..256");
    }
    for (const std::size_t index : group_indices) {
        if (index >= groups.size()) {
            throw std::invalid_argument(
                "projective bilinear-sum group index out of range");
        }
    }
    int worker_count = 1;
#ifdef NS_HAVE_OPENMP
    worker_count = std::min(
        threads,
        std::max(1, static_cast<int>(group_indices.size())));
#else
    static_cast<void>(threads);
#endif
    std::vector<SpectralIncrement> partials(
        static_cast<std::size_t>(worker_count),
        SpectralIncrement(state.waves.size()));
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
#ifdef NS_HAVE_OPENMP
#pragma omp parallel num_threads(worker_count) if(worker_count > 1)
    {
        const int worker = omp_get_thread_num();
        SpectralIncrement& partial = partials[
            static_cast<std::size_t>(worker)];
#pragma omp for schedule(dynamic, 1)
        for (std::ptrdiff_t position = 0;
             position < static_cast<std::ptrdiff_t>(group_indices.size());
             ++position) {
            const auto& group = groups[group_indices[
                static_cast<std::size_t>(position)]];
            for (const InteractionIndex interaction : group.interactions) {
                const auto [p, q, target] = interaction;
                const SpectralComplex coefficient = imaginary_unit *
                    wave_dot(state.waves[q], advecting[p]);
                for (std::size_t component = 0; component < 3;
                     ++component) {
                    partial[target][component] +=
                        coefficient * advected[q][component];
                }
            }
        }
    }
#else
    SpectralIncrement& partial = partials.front();
    for (const std::size_t index : group_indices) {
        for (const InteractionIndex interaction :
             groups[index].interactions) {
            const auto [p, q, target] = interaction;
            const SpectralComplex coefficient = imaginary_unit *
                wave_dot(state.waves[q], advecting[p]);
            for (std::size_t component = 0; component < 3;
                 ++component) {
                partial[target][component] +=
                    coefficient * advected[q][component];
            }
        }
    }
#endif
    SpectralIncrement result(state.waves.size());
    for (const SpectralIncrement& partial : partials) {
        for (std::size_t mode = 0; mode < result.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                result[mode][component] += partial[mode][component];
            }
        }
    }
    project(result, state);
    return result;
}

SpectralIncrement ProjectiveAdvectionDecomposition::vjp(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SpectralIncrement& output_cotangent) {
    ProjectiveBilinearCotangents cotangents = bilinear_vjp(
        state, group, state.velocity, state.velocity,
        output_cotangent);
    for (std::size_t mode = 0; mode < cotangents.advecting.size();
         ++mode) {
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            cotangents.advecting[mode][coordinate] +=
                cotangents.advected[mode][coordinate];
        }
    }
    return cotangents.advecting;
}

ProjectiveBilinearCotangents
ProjectiveAdvectionDecomposition::bilinear_vjp(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SpectralIncrement& advecting,
    const SpectralIncrement& advected,
    const SpectralIncrement& output_cotangent) {
    require_layout(state, advecting);
    require_layout(state, advected);
    require_layout(state, output_cotangent);
    SpectralIncrement cotangent = output_cotangent;
    project(cotangent, state);
    ProjectiveBilinearCotangents result{
        SpectralIncrement(state.waves.size()),
        SpectralIncrement(state.waves.size())};
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const ComplexVector& target_cotangent = cotangent[target];
        const SpectralComplex first_coefficient =
            minus_imaginary_unit *
            dot_hermitian(advected[q], target_cotangent);
        for (std::size_t component = 0; component < 3; ++component) {
            const SpectralReal wave_component =
                static_cast<SpectralReal>(
                    component == 0   ? state.waves[q].x
                    : component == 1 ? state.waves[q].y
                                     : state.waves[q].z);
            result.advecting[p][component] +=
                wave_component * first_coefficient;
        }
        const SpectralComplex second_coefficient =
                minus_imaginary_unit * std::conj(
                wave_dot(state.waves[q], advecting[p]));
        for (std::size_t component = 0; component < 3; ++component) {
            result.advected[q][component] +=
                second_coefficient * target_cotangent[component];
        }
    }
    project(result.advecting, state);
    project(result.advected, state);
    return result;
}

ProjectiveSquareFunctionMoment
ProjectiveAdvectionDecomposition::square_function(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    bool compute_gradient) {
    ProjectiveSquareFunctionMoment result;
    if (compute_gradient) {
        result.gradient.resize(state.waves.size());
    }
    SpectralIncrement component(state.waves.size());
    std::vector<std::size_t> touched_targets;
    touched_targets.reserve(state.waves.size());
    std::vector<std::size_t> target_generation(
        state.waves.size(), 0);
    std::size_t generation = 0;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const ProjectiveInteractionGroup& group : groups) {
        ++generation;
        touched_targets.clear();
        for (const InteractionIndex interaction : group.interactions) {
            const auto [p, q, target] = interaction;
            if (target_generation[target] != generation) {
                target_generation[target] = generation;
                component[target] = {};
                touched_targets.push_back(target);
            }
            const SpectralComplex coefficient = imaginary_unit *
                wave_dot(state.waves[q], state.velocity[p]);
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                component[target][coordinate] +=
                    coefficient * state.velocity[q][coordinate];
            }
        }
        for (const std::size_t target : touched_targets) {
            component[target] = project_divergence_free(
                state.waves[target], component[target]);
            result.norm2 += std::real(dot_hermitian(
                component[target], component[target]));
        }
        if (!compute_gradient) {
            continue;
        }
        for (const InteractionIndex interaction : group.interactions) {
            const auto [p, q, target] = interaction;
            ComplexVector target_cotangent = component[target];
            for (SpectralComplex& value : target_cotangent) {
                value *= 2.0L;
            }
            const SpectralComplex first_coefficient =
                minus_imaginary_unit *
                dot_hermitian(
                    state.velocity[q], target_cotangent);
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                const SpectralReal wave_component =
                    static_cast<SpectralReal>(
                        coordinate == 0   ? state.waves[q].x
                        : coordinate == 1 ? state.waves[q].y
                                          : state.waves[q].z);
                result.gradient[p][coordinate] +=
                    wave_component * first_coefficient;
            }
            const SpectralComplex second_coefficient =
                minus_imaginary_unit * std::conj(
                    wave_dot(state.waves[q], state.velocity[p]));
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                result.gradient[q][coordinate] +=
                    second_coefficient * target_cotangent[coordinate];
            }
        }
    }
    if (compute_gradient) {
        project(result.gradient, state);
    }
    return result;
}

ProjectiveSquareFunctionNorms
ProjectiveAdvectionDecomposition::square_function_norms(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    const std::vector<std::size_t>& group_indices,
    int threads) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective square-function workers must be 1..256");
    }
    for (const std::size_t index : group_indices) {
        if (index >= groups.size()) {
            throw std::invalid_argument(
                "projective square-function group index out of range");
        }
    }
    int worker_count = 1;
#ifdef NS_HAVE_OPENMP
    worker_count = std::min(
        threads,
        std::max(1, static_cast<int>(group_indices.size())));
#else
    static_cast<void>(threads);
#endif
    std::vector<ProjectiveSquareFunctionNorms> partials(
        static_cast<std::size_t>(worker_count));
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
#ifdef NS_HAVE_OPENMP
#pragma omp parallel num_threads(worker_count) if(worker_count > 1)
    {
        const int worker = omp_get_thread_num();
        auto& partial = partials[static_cast<std::size_t>(worker)];
        SpectralIncrement component(state.waves.size());
        std::vector<std::size_t> touched_targets;
        touched_targets.reserve(state.waves.size());
        std::vector<std::size_t> target_generation(
            state.waves.size(), 0);
        std::size_t generation = 0;
#pragma omp for schedule(dynamic, 1)
        for (std::ptrdiff_t position = 0;
             position < static_cast<std::ptrdiff_t>(group_indices.size());
             ++position) {
            ++generation;
            touched_targets.clear();
            const auto& group = groups[group_indices[
                static_cast<std::size_t>(position)]];
            for (const InteractionIndex interaction : group.interactions) {
                const auto [p, q, target] = interaction;
                if (target_generation[target] != generation) {
                    target_generation[target] = generation;
                    component[target] = {};
                    touched_targets.push_back(target);
                }
                const SpectralComplex coefficient = imaginary_unit *
                    wave_dot(state.waves[q], state.velocity[p]);
                for (std::size_t coordinate = 0; coordinate < 3;
                     ++coordinate) {
                    component[target][coordinate] +=
                        coefficient * state.velocity[q][coordinate];
                }
            }
            for (const std::size_t target : touched_targets) {
                component[target] = project_divergence_free(
                    state.waves[target], component[target]);
                const SpectralReal norm2 = std::real(dot_hermitian(
                    component[target], component[target]));
                const SpectralReal weight = static_cast<SpectralReal>(
                    norm_squared(state.waves[target]));
                partial.l2_norm2 += norm2;
                partial.h1_norm2 += weight * norm2;
                partial.h2_norm2 += weight * weight * norm2;
            }
        }
    }
#else
    auto& partial = partials.front();
    for (const std::size_t index : group_indices) {
        const SpectralIncrement component = evaluate(state, groups[index]);
        for (std::size_t target = 0; target < component.size(); ++target) {
            const SpectralReal norm2 = std::real(dot_hermitian(
                component[target], component[target]));
            const SpectralReal weight = static_cast<SpectralReal>(
                norm_squared(state.waves[target]));
            partial.l2_norm2 += norm2;
            partial.h1_norm2 += weight * norm2;
            partial.h2_norm2 += weight * weight * norm2;
        }
    }
#endif
    ProjectiveSquareFunctionNorms result;
    for (const auto& partial : partials) {
        result.l2_norm2 += partial.l2_norm2;
        result.h1_norm2 += partial.h1_norm2;
        result.h2_norm2 += partial.h2_norm2;
    }
    return result;
}

ProjectiveSquareFunctionMoment
ProjectiveAdvectionDecomposition::h1_square_function(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    const std::vector<std::size_t>& group_indices,
    bool compute_gradient) {
    for (const std::size_t index : group_indices) {
        if (index >= groups.size()) {
            throw std::invalid_argument(
                "projective H1 square-function group index out of range");
        }
    }
    ProjectiveSquareFunctionMoment result;
    if (compute_gradient) {
        result.gradient.resize(state.waves.size());
    }
    SpectralIncrement component(state.waves.size());
    std::vector<std::size_t> touched_targets;
    touched_targets.reserve(state.waves.size());
    std::vector<std::size_t> target_generation(
        state.waves.size(), 0);
    std::size_t generation = 0;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const std::size_t index : group_indices) {
        ++generation;
        touched_targets.clear();
        const ProjectiveInteractionGroup& group = groups[index];
        for (const InteractionIndex interaction : group.interactions) {
            const auto [p, q, target] = interaction;
            if (target_generation[target] != generation) {
                target_generation[target] = generation;
                component[target] = {};
                touched_targets.push_back(target);
            }
            const SpectralComplex coefficient = imaginary_unit *
                wave_dot(state.waves[q], state.velocity[p]);
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                component[target][coordinate] +=
                    coefficient * state.velocity[q][coordinate];
            }
        }
        for (const std::size_t target : touched_targets) {
            component[target] = project_divergence_free(
                state.waves[target], component[target]);
            const SpectralReal weight = static_cast<SpectralReal>(
                norm_squared(state.waves[target]));
            result.norm2 += weight * std::real(dot_hermitian(
                component[target], component[target]));
        }
        if (!compute_gradient) {
            continue;
        }
        for (const InteractionIndex interaction : group.interactions) {
            const auto [p, q, target] = interaction;
            ComplexVector target_cotangent = component[target];
            const SpectralReal cotangent_scale =
                2.0L * static_cast<SpectralReal>(
                    norm_squared(state.waves[target]));
            for (SpectralComplex& value : target_cotangent) {
                value *= cotangent_scale;
            }
            const SpectralComplex first_coefficient =
                minus_imaginary_unit * dot_hermitian(
                    state.velocity[q], target_cotangent);
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                const SpectralReal wave_component =
                    static_cast<SpectralReal>(
                        coordinate == 0   ? state.waves[q].x
                        : coordinate == 1 ? state.waves[q].y
                                          : state.waves[q].z);
                result.gradient[p][coordinate] +=
                    wave_component * first_coefficient;
            }
            const SpectralComplex second_coefficient =
                minus_imaginary_unit * std::conj(
                    wave_dot(state.waves[q], state.velocity[p]));
            for (std::size_t coordinate = 0; coordinate < 3;
                 ++coordinate) {
                result.gradient[q][coordinate] +=
                    second_coefficient * target_cotangent[coordinate];
            }
        }
    }
    if (compute_gradient) {
        project(result.gradient, state);
    }
    return result;
}

}  // namespace lemma
