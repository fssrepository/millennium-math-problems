#include "projective_quartic_diagonal_kernel.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <vector>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {
namespace {

using SparseField = std::map<std::size_t, ComplexVector>;

ComplexVector scaled(
    ComplexVector value,
    SpectralReal factor) {
    for (SpectralComplex& component : value) {
        component *= factor;
    }
    return value;
}

ComplexVector au_at(
    const SpectralState& state,
    std::size_t mode,
    int power = 1) {
    SpectralReal weight = 1.0L;
    const SpectralReal wave2 = static_cast<SpectralReal>(
        norm_squared(state.waves[mode]));
    for (int exponent = 0; exponent < power; ++exponent) {
        weight *= wave2;
    }
    return scaled(state.velocity[mode], weight);
}

ComplexVector sparse_value(
    const SparseField& field,
    std::size_t mode) {
    const auto position = field.find(mode);
    return position == field.end() ? ComplexVector{} : position->second;
}

void add_sparse(
    SparseField& target,
    std::size_t mode,
    const ComplexVector& value,
    SpectralReal factor = 1.0L) {
    ComplexVector& destination = target[mode];
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
        destination[coordinate] += factor * value[coordinate];
    }
}

void add_dense(
    SpectralIncrement& target,
    std::size_t mode,
    const ComplexVector& value,
    SpectralReal factor = 1.0L) {
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
        target[mode][coordinate] += factor * value[coordinate];
    }
}

SpectralReal sparse_pairing(
    const SparseField& left,
    const SparseField& right) {
    SpectralReal result = 0.0L;
    const SparseField* smaller = &left;
    const SparseField* larger = &right;
    if (left.size() > right.size()) {
        smaller = &right;
        larger = &left;
    }
    for (const auto& [mode, value] : *smaller) {
        const auto position = larger->find(mode);
        if (position != larger->end()) {
            result += std::real(dot_hermitian(
                value, position->second));
        }
    }
    return result;
}

SpectralReal sparse_pairing_au(
    const SpectralState& state,
    const SparseField& field,
    int au_power = 1) {
    SpectralReal result = 0.0L;
    for (const auto& [mode, value] : field) {
        result += std::real(dot_hermitian(
            au_at(state, mode, au_power), value));
    }
    return result;
}

SparseField laplacian_weight(
    const SpectralState& state,
    const SparseField& source) {
    SparseField result;
    for (const auto& [mode, value] : source) {
        result.emplace(
            mode,
            scaled(value, static_cast<SpectralReal>(
                norm_squared(state.waves[mode]))));
    }
    return result;
}

SparseField bilinear_full_full(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SpectralIncrement& advecting,
    const SpectralIncrement& advected) {
    SparseField result;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const SpectralComplex coefficient = imaginary_unit *
            wave_dot(state.waves[q], advecting[p]);
        ComplexVector& destination = result[target];
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            destination[coordinate] +=
                coefficient * advected[q][coordinate];
        }
    }
    for (auto& [mode, value] : result) {
        value = project_divergence_free(state.waves[mode], value);
    }
    return result;
}

SparseField bilinear_sparse_full(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SparseField& advecting,
    const SpectralIncrement& advected) {
    SparseField result;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const ComplexVector advecting_value = sparse_value(
            advecting, p);
        const SpectralComplex coefficient = imaginary_unit *
            wave_dot(state.waves[q], advecting_value);
        ComplexVector& destination = result[target];
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            destination[coordinate] +=
                coefficient * advected[q][coordinate];
        }
    }
    for (auto& [mode, value] : result) {
        value = project_divergence_free(state.waves[mode], value);
    }
    return result;
}

template <class AdvectingGetter,
          class AdvectedGetter,
          class AdvectingAdder,
          class AdvectedAdder>
void reverse_bilinear(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SparseField& output_cotangent,
    AdvectingGetter get_advecting,
    AdvectedGetter get_advected,
    AdvectingAdder add_advecting,
    AdvectedAdder add_advected) {
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const InteractionIndex interaction : group.interactions) {
        const auto [p, q, target] = interaction;
        const auto cotangent_position = output_cotangent.find(target);
        if (cotangent_position == output_cotangent.end()) {
            continue;
        }
        const ComplexVector target_cotangent =
            project_divergence_free(
                state.waves[target], cotangent_position->second);
        const ComplexVector advecting = get_advecting(p);
        const ComplexVector advected = get_advected(q);
        const SpectralComplex first_coefficient =
            minus_imaginary_unit *
            dot_hermitian(advected, target_cotangent);
        ComplexVector first{};
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            const SpectralReal wave_component =
                static_cast<SpectralReal>(
                    coordinate == 0   ? state.waves[q].x
                    : coordinate == 1 ? state.waves[q].y
                                      : state.waves[q].z);
            first[coordinate] =
                wave_component * first_coefficient;
        }
        add_advecting(p, first);
        const SpectralComplex second_coefficient =
            minus_imaginary_unit * std::conj(
                wave_dot(state.waves[q], advecting));
        ComplexVector second{};
        for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
            second[coordinate] =
                second_coefficient * target_cotangent[coordinate];
        }
        add_advected(q, second);
    }
}

struct SparseGroupMoment {
    SpectralReal bracket = 0.0L;
    SpectralReal d_enstrophy = 0.0L;
    SpectralReal d_palinstrophy = 0.0L;
};

SparseGroupMoment accumulate_group(
    const SpectralState& state,
    const ProjectiveInteractionGroup& group,
    const SpectralIncrement& au,
    SpectralReal enstrophy,
    SpectralReal palinstrophy,
    SpectralIncrement* gradient) {
    const SparseField b = bilinear_full_full(
        state, group, state.velocity, state.velocity);
    const SparseField ab = laplacian_weight(state, b);
    const SparseField transported = bilinear_full_full(
        state, group, state.velocity, au);
    const SparseField nested = bilinear_sparse_full(
        state, group, b, state.velocity);
    const SpectralReal stretching = sparse_pairing_au(state, b);
    const SpectralReal cross = sparse_pairing_au(state, ab);
    SparseGroupMoment moment;
    moment.bracket = -sparse_pairing(b, ab) +
        sparse_pairing(b, transported) -
        sparse_pairing_au(state, nested) +
        stretching * stretching / (2.0L * enstrophy) +
        3.0L * stretching * cross / (2.0L * palinstrophy);
    moment.d_enstrophy = -stretching * stretching /
        (2.0L * enstrophy * enstrophy);
    moment.d_palinstrophy = -3.0L * stretching * cross /
        (2.0L * palinstrophy * palinstrophy);
    if (gradient == nullptr) {
        return moment;
    }

    const SpectralReal d_stretching = stretching / enstrophy +
        3.0L * cross / (2.0L * palinstrophy);
    const SpectralReal d_cross =
        3.0L * stretching / (2.0L * palinstrophy);
    SparseField bar_b;
    SparseField bar_transported;
    SparseField bar_nested;
    SparseField bar_au;
    for (const auto& [mode, value] : b) {
        const ComplexVector ab_value = sparse_value(ab, mode);
        const ComplexVector transported_value = sparse_value(
            transported, mode);
        const ComplexVector au_value = au_at(state, mode);
        add_sparse(bar_b, mode, ab_value, -1.0L);
        add_sparse(bar_b, mode, transported_value);
        add_sparse(bar_b, mode, au_value, d_stretching);
        add_sparse(bar_transported, mode, value);
        add_sparse(bar_au, mode, value, d_stretching);

        ComplexVector bar_ab = scaled(value, -1.0L);
        add_sparse(
            bar_b, mode,
            scaled(au_value, d_cross * static_cast<SpectralReal>(
                norm_squared(state.waves[mode]))));
        add_sparse(
            bar_b, mode,
            scaled(bar_ab, static_cast<SpectralReal>(
                norm_squared(state.waves[mode]))));
        add_sparse(bar_au, mode, ab_value, d_cross);
    }
    for (const auto& [mode, value] : nested) {
        add_sparse(bar_au, mode, value, -1.0L);
        add_sparse(bar_nested, mode, au_at(state, mode), -1.0L);
    }

    reverse_bilinear(
        state, group, bar_nested,
        [&](std::size_t mode) { return sparse_value(b, mode); },
        [&](std::size_t mode) { return state.velocity[mode]; },
        [&](std::size_t mode, const ComplexVector& value) {
            add_sparse(bar_b, mode, value);
        },
        [&](std::size_t mode, const ComplexVector& value) {
            add_dense(*gradient, mode, value);
        });
    reverse_bilinear(
        state, group, bar_transported,
        [&](std::size_t mode) { return state.velocity[mode]; },
        [&](std::size_t mode) { return au[mode]; },
        [&](std::size_t mode, const ComplexVector& value) {
            add_dense(*gradient, mode, value);
        },
        [&](std::size_t mode, const ComplexVector& value) {
            add_sparse(bar_au, mode, value);
        });
    reverse_bilinear(
        state, group, bar_b,
        [&](std::size_t mode) { return state.velocity[mode]; },
        [&](std::size_t mode) { return state.velocity[mode]; },
        [&](std::size_t mode, const ComplexVector& value) {
            add_dense(*gradient, mode, value);
        },
        [&](std::size_t mode, const ComplexVector& value) {
            add_dense(*gradient, mode, value);
        });
    for (const auto& [mode, value] : bar_au) {
        add_dense(
            *gradient, mode, value,
            static_cast<SpectralReal>(norm_squared(state.waves[mode])));
    }
    return moment;
}

}  // namespace

ProjectiveQuarticDiagonalMoment
ProjectiveQuarticDiagonalKernel::evaluate(
    const SpectralState& state,
    const std::vector<ProjectiveInteractionGroup>& groups,
    SpectralReal enstrophy,
    SpectralReal palinstrophy,
    bool compute_gradient,
    int threads) {
    if (!(enstrophy > 0.0L) || !(palinstrophy > 0.0L)) {
        throw std::invalid_argument(
            "projective diagonal kernel requires positive Z and P");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective diagonal threads must be 1..256");
    }
    ProjectiveQuarticDiagonalMoment result;
    result.projective_shape_count = groups.size();
    SpectralIncrement au(state.waves.size());
    for (std::size_t mode = 0; mode < au.size(); ++mode) {
        au[mode] = au_at(state, mode);
    }
    int worker_count = 1;
#ifdef NS_HAVE_OPENMP
    worker_count = std::max(
        1, std::min(threads, static_cast<int>(groups.size())));
#endif
    std::vector<SpectralReal> partial_bracket(
        static_cast<std::size_t>(worker_count), 0.0L);
    std::vector<SpectralReal> partial_d_enstrophy(
        static_cast<std::size_t>(worker_count), 0.0L);
    std::vector<SpectralReal> partial_d_palinstrophy(
        static_cast<std::size_t>(worker_count), 0.0L);
    std::vector<SpectralIncrement> partial_gradients;
    if (compute_gradient) {
        partial_gradients.assign(
            static_cast<std::size_t>(worker_count),
            SpectralIncrement(state.waves.size()));
    }
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(worker_count) \
    if(worker_count > 1)
#endif
    for (std::ptrdiff_t group_index = 0;
         group_index < static_cast<std::ptrdiff_t>(groups.size());
         ++group_index) {
        int worker = 0;
#ifdef NS_HAVE_OPENMP
        worker = omp_get_thread_num();
#endif
        const std::size_t worker_index = static_cast<std::size_t>(worker);
        const SparseGroupMoment group_moment = accumulate_group(
            state, groups[static_cast<std::size_t>(group_index)], au,
            enstrophy, palinstrophy,
            compute_gradient
                ? &partial_gradients[worker_index] : nullptr);
        partial_bracket[worker_index] += group_moment.bracket;
        partial_d_enstrophy[worker_index] += group_moment.d_enstrophy;
        partial_d_palinstrophy[worker_index] +=
            group_moment.d_palinstrophy;
    }
    SpectralReal d_enstrophy = 0.0L;
    SpectralReal d_palinstrophy = 0.0L;
    if (compute_gradient) {
        result.gradient.resize(state.waves.size());
    }
    for (std::size_t worker = 0;
         worker < static_cast<std::size_t>(worker_count); ++worker) {
        result.bracket += partial_bracket[worker];
        d_enstrophy += partial_d_enstrophy[worker];
        d_palinstrophy += partial_d_palinstrophy[worker];
        if (compute_gradient) {
            for (std::size_t mode = 0;
                 mode < result.gradient.size(); ++mode) {
                add_dense(
                    result.gradient, mode,
                    partial_gradients[worker][mode]);
            }
        }
    }
    if (compute_gradient) {
        for (std::size_t mode = 0; mode < result.gradient.size(); ++mode) {
            const ComplexVector au = au_at(state, mode);
            const ComplexVector a2u = au_at(state, mode, 2);
            add_dense(result.gradient, mode, au, 2.0L * d_enstrophy);
            add_dense(
                result.gradient, mode, a2u,
                2.0L * d_palinstrophy);
            result.gradient[mode] = project_divergence_free(
                state.waves[mode], result.gradient[mode]);
        }
    }
    return result;
}

}  // namespace lemma
