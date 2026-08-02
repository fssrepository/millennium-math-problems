#include "local_quartic_closure_objective.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal pairing(const SpectralIncrement& left,
                     const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "local closure objective layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

void add_scaled(SpectralIncrement& target,
                const SpectralIncrement& source,
                SpectralReal scale) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "local closure gradient layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += scale * source[mode][component];
        }
    }
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source,
    int power = 1) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "local closure Laplacian layout mismatch");
    }
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        SpectralReal weight = 1.0L;
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (int exponent = 0; exponent < power; ++exponent) {
            weight *= wave2;
        }
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

struct ClosureGraph {
    SpectralIncrement au;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralIncrement transported_au;
    SpectralIncrement b_advects_u;
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal stretching = 0.0L;
    SpectralReal cross = 0.0L;
};

ClosureGraph build_graph(const SpectralDynamics& dynamics,
                         const SpectralState& state,
                         TriadSelection selection) {
    ClosureGraph graph;
    graph.au = laplacian_weight(state, state.velocity);
    graph.b = dynamics.advection_direct_partition(
        state, selection);
    graph.ab = laplacian_weight(state, graph.b);
    graph.transported_au =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, graph.au,
            selection);
    graph.b_advects_u =
        dynamics.advection_bilinear_direct_partition(
            state, graph.b, state.velocity,
            selection);
    graph.energy = pairing(state.velocity, state.velocity);
    graph.enstrophy = pairing(state.velocity, graph.au);
    graph.palinstrophy = pairing(graph.au, graph.au);
    graph.stretching = pairing(graph.au, graph.b);
    graph.cross = pairing(graph.ab, graph.au);
    return graph;
}

SpectralReal bracket_value(const ClosureGraph& graph) {
    if (!(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return 0.0L;
    }
    const SpectralReal negative_commutator = -pairing(
        graph.b, graph.ab) + pairing(
            graph.b, graph.transported_au);
    const SpectralReal advecting = -pairing(
        graph.au, graph.b_advects_u);
    return negative_commutator +
        graph.stretching * graph.stretching /
            (2.0L * graph.enstrophy) +
        3.0L * graph.stretching * graph.cross /
            (2.0L * graph.palinstrophy) +
        advecting;
}

SpectralIncrement bracket_gradient(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const ClosureGraph& graph,
    TriadSelection selection) {
    const std::size_t modes = state.waves.size();
    SpectralIncrement bar_u(modes);
    SpectralIncrement bar_au(modes);
    SpectralIncrement bar_b(modes);
    SpectralIncrement bar_ab(modes);
    SpectralIncrement bar_transported_au(modes);
    SpectralIncrement bar_b_advects_u(modes);

    add_scaled(bar_b, graph.ab, -1.0L);
    add_scaled(bar_ab, graph.b, -1.0L);
    add_scaled(bar_b, graph.transported_au, 1.0L);
    add_scaled(bar_transported_au, graph.b, 1.0L);
    add_scaled(bar_au, graph.b_advects_u, -1.0L);
    add_scaled(bar_b_advects_u, graph.au, -1.0L);

    const SpectralReal d_stretching =
        graph.stretching / graph.enstrophy +
        3.0L * graph.cross /
            (2.0L * graph.palinstrophy);
    const SpectralReal d_enstrophy =
        -graph.stretching * graph.stretching /
        (2.0L * graph.enstrophy * graph.enstrophy);
    const SpectralReal d_cross =
        3.0L * graph.stretching /
        (2.0L * graph.palinstrophy);
    const SpectralReal d_palinstrophy =
        -3.0L * graph.stretching * graph.cross /
        (2.0L * graph.palinstrophy * graph.palinstrophy);

    add_scaled(bar_au, graph.b, d_stretching);
    add_scaled(bar_b, graph.au, d_stretching);
    add_scaled(bar_u, graph.au, d_enstrophy);
    add_scaled(bar_au, state.velocity, d_enstrophy);
    add_scaled(bar_ab, graph.au, d_cross);
    add_scaled(bar_au, graph.ab, d_cross);
    add_scaled(bar_au, graph.au, 2.0L * d_palinstrophy);

    const BilinearAdvectionCotangents b_advects_u_cotangents =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, graph.b, state.velocity,
            bar_b_advects_u, selection);
    add_scaled(bar_b, b_advects_u_cotangents.advecting, 1.0L);
    add_scaled(bar_u, b_advects_u_cotangents.advected, 1.0L);

    const BilinearAdvectionCotangents transported_au_cotangents =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, state.velocity, graph.au,
            bar_transported_au, selection);
    add_scaled(bar_u, transported_au_cotangents.advecting, 1.0L);
    add_scaled(bar_au, transported_au_cotangents.advected, 1.0L);

    add_scaled(bar_b, laplacian_weight(state, bar_ab), 1.0L);
    add_scaled(
        bar_u,
        dynamics.advection_vjp_direct_partition(
            state, bar_b, selection),
        1.0L);
    add_scaled(bar_u, laplacian_weight(state, bar_au), 1.0L);
    return bar_u;
}

SpectralIncrement stretching_gradient(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const ClosureGraph& graph,
    TriadSelection selection) {
    SpectralIncrement result = graph.ab;
    add_scaled(
        result,
        dynamics.advection_vjp_direct_partition(
            state, graph.au, selection),
        1.0L);
    return result;
}

}  // namespace

LocalQuarticClosureObjective::LocalQuarticClosureObjective(
    const SpectralDynamics& dynamics, TriadSelection selection)
    : dynamics_(dynamics), selection_(selection) {}

LocalQuarticClosureObjectiveValue LocalQuarticClosureObjective::evaluate(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    LocalQuarticClosureObjectiveValue result;
    result.energy = graph.energy;
    result.enstrophy = graph.enstrophy;
    result.palinstrophy = graph.palinstrophy;
    result.signed_stretching = graph.stretching;
    result.palinstrophy_cross = graph.cross;
    result.negative_commutator_pairing = -pairing(
        graph.b, graph.ab) + pairing(
            graph.b, graph.transported_au);
    result.advecting_slot = -pairing(
        graph.au, graph.b_advects_u);
    result.signed_two_entry_bracket = bracket_value(graph);
    if (result.energy > 0.0L && result.enstrophy > 0.0L &&
        result.palinstrophy > 0.0L) {
        result.candidate_scale =
            std::pow(result.enstrophy, 1.75L) *
            result.palinstrophy /
            std::pow(result.energy, 0.25L);
        result.constant_ratio = std::abs(
            result.signed_two_entry_bracket) /
            result.candidate_scale;
        result.signed_constant_ratio =
            result.signed_two_entry_bracket /
            result.candidate_scale;
        result.squared_constant_ratio =
            result.constant_ratio * result.constant_ratio;
        result.lqc3_target_scale =
            std::pow(result.enstrophy, 1.25L) *
            std::pow(result.palinstrophy, 0.75L);
        result.lqc3_target_ratio = std::abs(
            result.signed_two_entry_bracket) /
            result.lqc3_target_scale;
        result.signed_lqc3_target_ratio =
            result.signed_two_entry_bracket /
            result.lqc3_target_scale;
        result.squared_lqc3_target_ratio =
            result.lqc3_target_ratio * result.lqc3_target_ratio;
        result.initial_frequency = std::sqrt(
            result.enstrophy / result.energy);
        result.initial_ep_shift = result.energy * result.palinstrophy;
        const SpectralReal stretching2 =
            result.signed_stretching * result.signed_stretching;
        const SpectralReal stretching4 = stretching2 * stretching2;
        result.local_polynomial_numerator =
            4.0L * result.signed_stretching * stretching2 *
            result.enstrophy * result.palinstrophy *
            result.signed_two_entry_bracket;
        result.local_polynomial_denominator =
            result.initial_frequency *
            (stretching4 * result.enstrophy * result.enstrophy *
                 result.palinstrophy +
             result.initial_ep_shift *
                 result.enstrophy * result.enstrophy *
                 result.enstrophy *
                 result.palinstrophy * result.palinstrophy *
                 result.palinstrophy * result.palinstrophy);
        if (result.local_polynomial_denominator > 0.0L) {
            result.signed_local_sld_ratio =
                result.local_polynomial_numerator /
                result.local_polynomial_denominator;
        }
        result.normalized_stretching_ratio =
            result.signed_stretching /
            (std::pow(result.energy, 0.25L) *
             std::pow(result.enstrophy, 0.25L) *
             result.palinstrophy);
        const SpectralReal normalized2 =
            result.normalized_stretching_ratio *
            result.normalized_stretching_ratio;
        const SpectralReal normalized3 = normalized2 *
            result.normalized_stretching_ratio;
        const SpectralReal normalized4 = normalized2 * normalized2;
        result.signed_shape_factor =
            4.0L * normalized3 / (1.0L + normalized4);
        result.factorized_local_sld_ratio =
            result.signed_constant_ratio * result.signed_shape_factor;
        const SpectralReal factorization_scale = std::max(
            {std::abs(result.signed_local_sld_ratio),
             std::abs(result.factorized_local_sld_ratio), 1e-30L});
        result.factorization_relative_error = std::abs(
            result.signed_local_sld_ratio -
            result.factorized_local_sld_ratio) /
            factorization_scale;
    }
    result.finite = std::isfinite(result.squared_constant_ratio) &&
        std::isfinite(result.squared_lqc3_target_ratio) &&
        std::isfinite(result.signed_lqc3_target_ratio) &&
        std::isfinite(result.signed_local_sld_ratio) &&
        std::isfinite(result.factorization_relative_error) &&
        result.candidate_scale > 0.0L;
    return result;
}

SpectralIncrement
LocalQuarticClosureObjective::signed_local_sld_ratio_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    const SpectralReal e = graph.energy;
    const SpectralReal z = graph.enstrophy;
    const SpectralReal p = graph.palinstrophy;
    const SpectralReal s = graph.stretching;
    const SpectralReal f = bracket_value(graph);
    if (!(e > 0.0L) || !(z > 0.0L) || !(p > 0.0L)) {
        return result;
    }
    const SpectralReal s2 = s * s;
    const SpectralReal s3 = s2 * s;
    const SpectralReal s4 = s2 * s2;
    const SpectralReal p2 = p * p;
    const SpectralReal p4 = p2 * p2;
    const SpectralReal first = s4 * z * z * p;
    const SpectralReal second = e * z * z * z * p * p4;
    const SpectralReal sum = first + second;
    const SpectralReal frequency = std::sqrt(z / e);
    const SpectralReal denominator = frequency * sum;
    if (!(denominator > 0.0L) || !std::isfinite(denominator)) {
        return result;
    }
    const SpectralReal numerator = 4.0L * s3 * z * p * f;
    const SpectralReal denominator2 = denominator * denominator;
    const SpectralReal numerator_f = 4.0L * s3 * z * p;
    const SpectralReal numerator_s = 12.0L * s2 * z * p * f;
    const SpectralReal numerator_z = 4.0L * s3 * p * f;
    const SpectralReal numerator_p = 4.0L * s3 * z * f;
    const SpectralReal denominator_e = frequency *
        (z * z * z * p * p4 - sum / (2.0L * e));
    const SpectralReal denominator_z = frequency *
        (2.0L * s4 * z * p + 3.0L * e * z * z * p * p4 +
         sum / (2.0L * z));
    const SpectralReal denominator_p = frequency *
        (s4 * z * z + 5.0L * e * z * z * z * p4);
    const SpectralReal denominator_s = frequency *
        (4.0L * s3 * z * z * p);
    auto quotient_partial = [&](SpectralReal numerator_partial,
                                SpectralReal denominator_partial) {
        return (numerator_partial * denominator -
                numerator * denominator_partial) / denominator2;
    };
    const SpectralReal weight_f = numerator_f / denominator;
    const SpectralReal weight_s = quotient_partial(
        numerator_s, denominator_s);
    const SpectralReal weight_e = quotient_partial(
        0.0L, denominator_e);
    const SpectralReal weight_z = quotient_partial(
        numerator_z, denominator_z);
    const SpectralReal weight_p = quotient_partial(
        numerator_p, denominator_p);
    result = bracket_gradient(dynamics_, state, graph, selection_);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= weight_f;
        }
    }
    add_scaled(
        result, stretching_gradient(
            dynamics_, state, graph, selection_), weight_s);
    add_scaled(result, state.velocity, 2.0L * weight_e);
    add_scaled(result, graph.au, 2.0L * weight_z);
    add_scaled(
        result, laplacian_weight(state, graph.au), 2.0L * weight_p);
    return result;
}

SpectralReal LocalQuarticClosureObjective::frozen_signed_local_sld_ratio(
    const SpectralState& state,
    SpectralReal initial_frequency,
    SpectralReal initial_ep_shift) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    const SpectralReal z = graph.enstrophy;
    const SpectralReal p = graph.palinstrophy;
    const SpectralReal s = graph.stretching;
    if (!(initial_frequency > 0.0L) ||
        !(initial_ep_shift > 0.0L) || !(z > 0.0L) || !(p > 0.0L)) {
        return 0.0L;
    }
    const SpectralReal s2 = s * s;
    const SpectralReal s3 = s2 * s;
    const SpectralReal s4 = s2 * s2;
    const SpectralReal p2 = p * p;
    const SpectralReal p4 = p2 * p2;
    const SpectralReal denominator = initial_frequency *
        (s4 * z * z * p +
         initial_ep_shift * z * z * z * p4);
    return denominator > 0.0L
        ? 4.0L * s3 * z * p * bracket_value(graph) / denominator
        : 0.0L;
}

SpectralIncrement
LocalQuarticClosureObjective::frozen_signed_local_sld_ratio_gradient(
    const SpectralState& state,
    SpectralReal initial_frequency,
    SpectralReal initial_ep_shift) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    const SpectralReal z = graph.enstrophy;
    const SpectralReal p = graph.palinstrophy;
    const SpectralReal s = graph.stretching;
    const SpectralReal f = bracket_value(graph);
    if (!(initial_frequency > 0.0L) ||
        !(initial_ep_shift > 0.0L) || !(z > 0.0L) || !(p > 0.0L)) {
        return result;
    }
    const SpectralReal s2 = s * s;
    const SpectralReal s3 = s2 * s;
    const SpectralReal s4 = s2 * s2;
    const SpectralReal p2 = p * p;
    const SpectralReal p3 = p2 * p;
    const SpectralReal p4 = p2 * p2;
    const SpectralReal inner = s4 * z * z * p +
        initial_ep_shift * z * z * z * p4;
    const SpectralReal denominator = initial_frequency * inner;
    if (!(denominator > 0.0L) || !std::isfinite(denominator)) {
        return result;
    }
    const SpectralReal numerator = 4.0L * s3 * z * p * f;
    const SpectralReal denominator2 = denominator * denominator;
    const SpectralReal numerator_f = 4.0L * s3 * z * p;
    const SpectralReal numerator_s = 12.0L * s2 * z * p * f;
    const SpectralReal numerator_z = 4.0L * s3 * p * f;
    const SpectralReal numerator_p = 4.0L * s3 * z * f;
    const SpectralReal denominator_s = initial_frequency *
        (4.0L * s3 * z * z * p);
    const SpectralReal denominator_z = initial_frequency *
        (2.0L * s4 * z * p +
         3.0L * initial_ep_shift * z * z * p4);
    const SpectralReal denominator_p = initial_frequency *
        (s4 * z * z +
         4.0L * initial_ep_shift * z * z * z * p3);
    auto quotient_partial = [&](SpectralReal numerator_partial,
                                SpectralReal denominator_partial) {
        return (numerator_partial * denominator -
                numerator * denominator_partial) / denominator2;
    };
    const SpectralReal weight_f = numerator_f / denominator;
    const SpectralReal weight_s = quotient_partial(
        numerator_s, denominator_s);
    const SpectralReal weight_z = quotient_partial(
        numerator_z, denominator_z);
    const SpectralReal weight_p = quotient_partial(
        numerator_p, denominator_p);
    result = bracket_gradient(dynamics_, state, graph, selection_);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= weight_f;
        }
    }
    add_scaled(
        result, stretching_gradient(
            dynamics_, state, graph, selection_), weight_s);
    add_scaled(result, graph.au, 2.0L * weight_z);
    add_scaled(
        result, laplacian_weight(state, graph.au), 2.0L * weight_p);
    return result;
}

SpectralIncrement
LocalQuarticClosureObjective::two_entry_bracket_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    if (!(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return SpectralIncrement(state.waves.size());
    }
    return bracket_gradient(dynamics_, state, graph, selection_);
}

SpectralIncrement
LocalQuarticClosureObjective::signed_stretching_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    if (!(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return SpectralIncrement(state.waves.size());
    }
    return stretching_gradient(
        dynamics_, state, graph, selection_);
}

SpectralIncrement
LocalQuarticClosureObjective::palinstrophy_cross_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    if (!(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return SpectralIncrement(state.waves.size());
    }
    SpectralIncrement result = laplacian_weight(state, graph.ab);
    add_scaled(
        result,
        dynamics_.advection_vjp_direct_partition(
            state, laplacian_weight(state, graph.au), selection_),
        1.0L);
    return result;
}

SpectralIncrement
LocalQuarticClosureObjective::squared_constant_ratio_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    if (!(graph.energy > 0.0L) || !(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal bracket = bracket_value(graph);
    const SpectralReal factor = std::sqrt(graph.energy) /
        (std::pow(graph.enstrophy, 3.5L) *
         graph.palinstrophy * graph.palinstrophy);
    const SpectralReal objective = bracket * bracket * factor;
    result = bracket_gradient(dynamics_, state, graph, selection_);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= 2.0L * bracket * factor;
        }
    }
    add_scaled(result, state.velocity, objective / graph.energy);
    add_scaled(result, graph.au, -7.0L * objective / graph.enstrophy);
    add_scaled(
        result, laplacian_weight(state, graph.au),
        -4.0L * objective / graph.palinstrophy);
    return result;
}

SpectralIncrement
LocalQuarticClosureObjective::signed_constant_ratio_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    if (!(graph.energy > 0.0L) || !(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal bracket = bracket_value(graph);
    const SpectralReal factor = std::pow(graph.energy, 0.25L) /
        (std::pow(graph.enstrophy, 1.75L) * graph.palinstrophy);
    const SpectralReal objective = bracket * factor;
    result = bracket_gradient(
        dynamics_, state, graph, selection_);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= factor;
        }
    }
    add_scaled(
        result, state.velocity, 0.5L * objective / graph.energy);
    add_scaled(
        result, graph.au, -3.5L * objective / graph.enstrophy);
    add_scaled(
        result, laplacian_weight(state, graph.au),
        -2.0L * objective / graph.palinstrophy);
    return result;
}

SpectralIncrement
LocalQuarticClosureObjective::squared_lqc3_target_ratio_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    if (!(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal bracket = bracket_value(graph);
    const SpectralReal factor = 1.0L /
        (std::pow(graph.enstrophy, 2.5L) *
         std::pow(graph.palinstrophy, 1.5L));
    const SpectralReal objective = bracket * bracket * factor;
    result = bracket_gradient(dynamics_, state, graph, selection_);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= 2.0L * bracket * factor;
        }
    }
    add_scaled(result, graph.au,
               -5.0L * objective / graph.enstrophy);
    add_scaled(result, laplacian_weight(state, graph.au),
               -3.0L * objective / graph.palinstrophy);
    return result;
}

SpectralIncrement
LocalQuarticClosureObjective::signed_lqc3_target_ratio_gradient(
    const SpectralState& state) const {
    const ClosureGraph graph = build_graph(dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    if (!(graph.enstrophy > 0.0L) ||
        !(graph.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal bracket = bracket_value(graph);
    const SpectralReal factor = 1.0L /
        (std::pow(graph.enstrophy, 1.25L) *
         std::pow(graph.palinstrophy, 0.75L));
    const SpectralReal objective = bracket * factor;
    result = bracket_gradient(dynamics_, state, graph, selection_);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= factor;
        }
    }
    add_scaled(result, graph.au,
               -2.5L * objective / graph.enstrophy);
    add_scaled(result, laplacian_weight(state, graph.au),
               -1.5L * objective / graph.palinstrophy);
    return result;
}

}  // namespace lemma
