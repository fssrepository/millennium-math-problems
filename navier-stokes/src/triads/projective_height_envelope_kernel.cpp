#include "projective_height_envelope_kernel.hpp"

#include "projective_height_shell_partition.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

SpectralReal pairing(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal sign(SpectralReal value) {
    return value > 0.0L ? 1.0L : (value < 0.0L ? -1.0L : 0.0L);
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

struct ShellGraph {
    const ProjectiveInteractionGroup* group = nullptr;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralIncrement transported_au;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
};

struct PairGraph {
    std::array<SpectralReal, 5> components{};
    std::array<SpectralReal, 5> signs{};
};

void accumulate_component(
    ProjectiveHeightEnvelopeMoment& result,
    PairGraph& pair,
    std::size_t component,
    SpectralReal value) {
    pair.components[component] = value;
    pair.signs[component] = sign(value);
    result.absolute_component_sums[component] += std::abs(value);
    result.absolute_component_envelope += std::abs(value);
}

void accumulate_dynamic_components(
    ProjectiveHeightEnvelopeMoment& result,
    PairGraph& pair,
    SpectralReal outer,
    SpectralReal advected,
    SpectralReal nested,
    bool pair_outer_and_advected,
    bool pair_nested_with_commutator) {
    pair.components[0] = outer;
    pair.components[1] = advected;
    pair.components[2] = nested;
    if (pair_nested_with_commutator) {
        const SpectralReal dynamic_sign = sign(
            outer + advected + nested);
        pair.signs[0] = dynamic_sign;
        pair.signs[1] = dynamic_sign;
        pair.signs[2] = dynamic_sign;
        result.absolute_component_sums[0] +=
            std::abs(outer + advected + nested);
        result.absolute_component_envelope +=
            std::abs(outer + advected + nested);
        return;
    }
    pair.signs[2] = sign(nested);
    result.absolute_component_sums[2] += std::abs(nested);
    result.absolute_component_envelope += std::abs(nested);
    if (pair_outer_and_advected) {
        const SpectralReal paired_sign = sign(outer + advected);
        pair.signs[0] = paired_sign;
        pair.signs[1] = paired_sign;
        result.absolute_component_sums[0] +=
            std::abs(outer + advected);
        result.absolute_component_envelope +=
            std::abs(outer + advected);
        return;
    }
    pair.signs[0] = sign(outer);
    pair.signs[1] = sign(advected);
    result.absolute_component_sums[0] += std::abs(outer);
    result.absolute_component_sums[1] += std::abs(advected);
    result.absolute_component_envelope +=
        std::abs(outer) + std::abs(advected);
}

}  // namespace

ProjectiveHeightEnvelopeMoment ProjectiveHeightEnvelopeKernel::evaluate(
    const SpectralState& state,
    TriadSelection selection,
    SpectralReal enstrophy,
    SpectralReal palinstrophy,
    bool compute_gradient,
    int threads,
    bool pair_outer_and_advected_commutator,
    bool pair_nested_with_commutator) {
    if (!(enstrophy > 0.0L) || !(palinstrophy > 0.0L)) {
        throw std::invalid_argument(
            "projective height envelope requires positive Z and P");
    }
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective height envelope threads must be 1..256");
    }
    const auto& groups = ProjectiveAdvectionDecomposition::group(
        state, selection);
    const auto partition = ProjectiveHeightShellPartition::build(groups);
    ProjectiveHeightEnvelopeMoment result;
    result.height_shell_count = partition.size();
    result.height_pair_count =
        partition.size() * (partition.size() + 1) / 2;
    if (partition.empty()) {
        if (compute_gradient) {
            result.gradient.resize(state.waves.size());
            result.outer_gradient.resize(state.waves.size());
        }
        return result;
    }

    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const SpectralIncrement a2u = laplacian_weight(state, au);
    std::vector<ShellGraph> shells(partition.size());
    for (std::size_t index = 0; index < partition.size(); ++index) {
        ShellGraph& shell = shells[index];
        const auto& aggregate =
            ProjectiveAdvectionDecomposition::aggregate_family(
                state, selection, partition[index].primitive_shapes);
        shell.group = &aggregate.front();
        shell.b = ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
            state, groups, partition[index].group_indices,
            state.velocity, state.velocity, threads);
        shell.ab = laplacian_weight(state, shell.b);
        result.outer_h1_sum += pairing(shell.b, shell.ab);
        shell.transported_au =
            ProjectiveAdvectionDecomposition::evaluate_bilinear_sum(
                state, groups, partition[index].group_indices,
                state.velocity, au, threads);
        shell.stretching = pairing(au, shell.b);
        shell.palinstrophy_cross = pairing(shell.ab, au);
    }

    const std::size_t shell_count = shells.size();
    std::vector<PairGraph> pairs(shell_count * shell_count);
    std::vector<SpectralIncrement> nested;
    if (compute_gradient) {
        nested.assign(
            shell_count * shell_count,
            SpectralIncrement(state.waves.size()));
    }
    auto nested_value = [&](std::size_t advector, std::size_t source) {
        SpectralIncrement value =
            ProjectiveAdvectionDecomposition::evaluate_bilinear(
                state, *shells[advector].group,
                shells[source].b, state.velocity, threads);
        if (compute_gradient) {
            nested[advector * shell_count + source] = value;
        }
        return pairing(au, value);
    };

    for (std::size_t first = 0; first < shell_count; ++first) {
        for (std::size_t second = first; second < shell_count; ++second) {
            PairGraph& pair = pairs[first * shell_count + second];
            const ShellGraph& left = shells[first];
            const ShellGraph& right = shells[second];
            if (first == second) {
                accumulate_dynamic_components(
                    result, pair,
                    -pairing(left.b, left.ab),
                    pairing(left.b, left.transported_au),
                    -nested_value(first, first),
                    pair_outer_and_advected_commutator,
                    pair_nested_with_commutator);
                accumulate_component(
                    result, pair, 3,
                    left.stretching * left.stretching /
                        (2.0L * enstrophy));
                accumulate_component(
                    result, pair, 4,
                    3.0L * left.stretching *
                        left.palinstrophy_cross /
                        (2.0L * palinstrophy));
            } else {
                accumulate_dynamic_components(
                    result, pair,
                    -pairing(left.b, right.ab) -
                        pairing(right.b, left.ab),
                    pairing(left.b, right.transported_au) +
                        pairing(right.b, left.transported_au),
                    -nested_value(first, second) -
                        nested_value(second, first),
                    pair_outer_and_advected_commutator,
                    pair_nested_with_commutator);
                accumulate_component(
                    result, pair, 3,
                    left.stretching * right.stretching /
                        enstrophy);
                accumulate_component(
                    result, pair, 4,
                    3.0L *
                        (left.stretching *
                             right.palinstrophy_cross +
                         right.stretching *
                             left.palinstrophy_cross) /
                        (2.0L * palinstrophy));
            }
        }
    }
    if (!compute_gradient) {
        return result;
    }

    result.gradient.resize(state.waves.size());
    result.outer_gradient.resize(state.waves.size());
    SpectralIncrement bar_au(state.waves.size());
    std::vector<SpectralIncrement> bar_b(
        shell_count, SpectralIncrement(state.waves.size()));
    std::vector<SpectralIncrement> bar_transported(
        shell_count, SpectralIncrement(state.waves.size()));
    std::vector<SpectralReal> bar_stretching(shell_count, 0.0L);
    std::vector<SpectralReal> bar_cross(shell_count, 0.0L);
    SpectralReal d_enstrophy = 0.0L;
    SpectralReal d_palinstrophy = 0.0L;

    auto reverse_nested = [&](std::size_t advector,
                              std::size_t source,
                              SpectralReal component_sign) {
        const SpectralReal coefficient = -component_sign;
        const SpectralIncrement& value =
            nested[advector * shell_count + source];
        add_scaled(bar_au, value, coefficient);
        const ProjectiveBilinearCotangents cotangents =
            ProjectiveAdvectionDecomposition::bilinear_vjp(
                state, *shells[advector].group,
                shells[source].b, state.velocity, au, threads);
        add_scaled(
            bar_b[source], cotangents.advecting, coefficient);
        add_scaled(
            result.gradient, cotangents.advected, coefficient);
    };

    for (std::size_t first = 0; first < shell_count; ++first) {
        for (std::size_t second = first; second < shell_count; ++second) {
            const PairGraph& pair = pairs[first * shell_count + second];
            const SpectralReal outer_sign = pair.signs[0];
            const SpectralReal transported_sign = pair.signs[1];
            const SpectralReal nested_sign = pair.signs[2];
            const SpectralReal enstrophy_sign = pair.signs[3];
            const SpectralReal palinstrophy_sign = pair.signs[4];
            if (first == second) {
                add_scaled(
                    bar_b[first], shells[first].ab,
                    -2.0L * outer_sign);
                add_scaled(
                    bar_b[first], shells[first].transported_au,
                    transported_sign);
                add_scaled(
                    bar_transported[first], shells[first].b,
                    transported_sign);
                reverse_nested(first, first, nested_sign);
                bar_stretching[first] += enstrophy_sign *
                    shells[first].stretching / enstrophy;
                d_enstrophy -= enstrophy_sign *
                    shells[first].stretching *
                    shells[first].stretching /
                    (2.0L * enstrophy * enstrophy);
                bar_stretching[first] += palinstrophy_sign *
                    3.0L * shells[first].palinstrophy_cross /
                    (2.0L * palinstrophy);
                bar_cross[first] += palinstrophy_sign *
                    3.0L * shells[first].stretching /
                    (2.0L * palinstrophy);
                d_palinstrophy -= palinstrophy_sign *
                    pair.components[4] / palinstrophy;
            } else {
                add_scaled(
                    bar_b[first], shells[second].ab,
                    -2.0L * outer_sign);
                add_scaled(
                    bar_b[second], shells[first].ab,
                    -2.0L * outer_sign);
                add_scaled(
                    bar_b[first], shells[second].transported_au,
                    transported_sign);
                add_scaled(
                    bar_b[second], shells[first].transported_au,
                    transported_sign);
                add_scaled(
                    bar_transported[first], shells[second].b,
                    transported_sign);
                add_scaled(
                    bar_transported[second], shells[first].b,
                    transported_sign);
                reverse_nested(first, second, nested_sign);
                reverse_nested(second, first, nested_sign);
                bar_stretching[first] += enstrophy_sign *
                    shells[second].stretching / enstrophy;
                bar_stretching[second] += enstrophy_sign *
                    shells[first].stretching / enstrophy;
                d_enstrophy -= enstrophy_sign *
                    shells[first].stretching *
                    shells[second].stretching /
                    (enstrophy * enstrophy);
                bar_stretching[first] += palinstrophy_sign *
                    3.0L * shells[second].palinstrophy_cross /
                    (2.0L * palinstrophy);
                bar_stretching[second] += palinstrophy_sign *
                    3.0L * shells[first].palinstrophy_cross /
                    (2.0L * palinstrophy);
                bar_cross[first] += palinstrophy_sign *
                    3.0L * shells[second].stretching /
                    (2.0L * palinstrophy);
                bar_cross[second] += palinstrophy_sign *
                    3.0L * shells[first].stretching /
                    (2.0L * palinstrophy);
                d_palinstrophy -= palinstrophy_sign *
                    pair.components[4] / palinstrophy;
            }
        }
    }

    for (std::size_t shell = 0; shell < shell_count; ++shell) {
        add_scaled(
            bar_au, shells[shell].b, bar_stretching[shell]);
        add_scaled(
            bar_b[shell], au, bar_stretching[shell]);
        add_scaled(
            bar_b[shell], a2u, bar_cross[shell]);
        add_scaled(
            bar_au, shells[shell].ab, bar_cross[shell]);

        const ProjectiveBilinearCotangents transported_cotangents =
            ProjectiveAdvectionDecomposition::bilinear_vjp(
                state, *shells[shell].group,
                state.velocity, au, bar_transported[shell], threads);
        add_scaled(
            result.gradient, transported_cotangents.advecting, 1.0L);
        add_scaled(
            bar_au, transported_cotangents.advected, 1.0L);
        add_scaled(
            result.gradient,
            ProjectiveAdvectionDecomposition::vjp(
                state, *shells[shell].group, bar_b[shell], threads),
            1.0L);
        add_scaled(
            result.outer_gradient,
            ProjectiveAdvectionDecomposition::vjp(
                state, *shells[shell].group,
                shells[shell].ab, threads),
            2.0L);
    }
    add_scaled(
        result.gradient,
        laplacian_weight(state, bar_au), 1.0L);
    add_scaled(result.gradient, au, 2.0L * d_enstrophy);
    add_scaled(result.gradient, a2u, 2.0L * d_palinstrophy);
    return result;
}

}  // namespace lemma
