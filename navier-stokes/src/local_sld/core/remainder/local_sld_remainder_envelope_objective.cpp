#include "local_sld_remainder_envelope_objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal pairing(const SpectralIncrement& left,
                     const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "remainder envelope layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source,
    int power) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "remainder envelope Laplacian mismatch");
    }
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        SpectralReal weight = 1.0L;
        if (power >= 0) {
            for (int index = 0; index < power; ++index) {
                weight *= wave2;
            }
        } else {
            for (int index = 0; index < -power; ++index) {
                weight /= wave2;
            }
        }
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

void add_scaled(SpectralIncrement& target,
                const SpectralIncrement& source,
                SpectralReal scale) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "remainder envelope gradient mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += scale * source[mode][component];
        }
    }
}

SpectralIncrement difference(const SpectralIncrement& left,
                             const SpectralIncrement& right) {
    SpectralIncrement result = left;
    add_scaled(result, right, -1.0L);
    return result;
}

struct EnvelopeGraph {
    SpectralIncrement au;
    SpectralIncrement aau;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralIncrement d;
    SpectralIncrement inverse_d;
    SpectralReal z = 0.0L;
    SpectralReal p = 0.0L;
    SpectralReal s = 0.0L;
    SpectralReal h3 = 0.0L;
    SpectralReal j = 0.0L;
    SpectralReal hminus1 = 0.0L;
    SpectralReal c = 0.0L;
};

EnvelopeGraph build_graph(const SpectralDynamics& dynamics,
                          const SpectralState& state,
                          TriadSelection selection) {
    EnvelopeGraph graph;
    graph.au = laplacian_weight(state, state.velocity, 1);
    graph.aau = laplacian_weight(state, graph.au, 1);
    graph.b = dynamics.advection_direct_partition(state, selection);
    graph.ab = laplacian_weight(state, graph.b, 1);
    const SpectralIncrement transported_au =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, graph.au, selection);
    const BilinearAdvectionCotangents nested_cotangents =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, state.velocity, state.velocity,
            graph.au, selection);
    graph.d = difference(
        transported_au, nested_cotangents.advecting);
    graph.inverse_d = laplacian_weight(state, graph.d, -1);
    graph.z = pairing(state.velocity, graph.au);
    graph.p = pairing(graph.au, graph.au);
    graph.s = pairing(graph.au, graph.b);
    graph.h3 = pairing(graph.au, graph.aau);
    graph.j = pairing(graph.au, graph.d);
    graph.hminus1 = pairing(graph.d, graph.inverse_d);
    graph.c = 3.0L * graph.s / (4.0L * graph.p);
    return graph;
}

SpectralIncrement d_vjp(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const EnvelopeGraph& graph,
    const SpectralIncrement& cotangent,
    TriadSelection selection) {
    const BilinearAdvectionCotangents transported =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, state.velocity, graph.au,
            cotangent, selection);
    const SpectralIncrement cotangent_advects_u =
        dynamics.advection_bilinear_direct_partition(
            state, cotangent, state.velocity, selection);
    const BilinearAdvectionCotangents nested =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, cotangent, state.velocity,
            graph.au, selection);
    SpectralIncrement result = transported.advecting;
    add_scaled(result, nested.advected, -1.0L);
    SpectralIncrement au_cotangent = transported.advected;
    add_scaled(au_cotangent, cotangent_advects_u, -1.0L);
    add_scaled(result, laplacian_weight(state, au_cotangent, 1), 1.0L);
    return result;
}

LocalSldRemainderEnvelopeValue value_from_graph(
    const EnvelopeGraph& graph) {
    LocalSldRemainderEnvelopeValue value;
    value.enstrophy = graph.z;
    value.palinstrophy = graph.p;
    value.stretching = graph.s;
    value.projection_coefficient = graph.c;
    value.enstrophy_normalization = graph.s * graph.s /
        (2.0L * graph.z);
    value.projected_h3_correction = graph.c * graph.c * graph.h3;
    value.projected_commutator_pairing = graph.c * graph.j;
    value.commutator_hminus1_norm2 = graph.hminus1;
    value.upper_envelope = value.enstrophy_normalization +
        value.projected_h3_correction +
        value.projected_commutator_pairing +
        0.25L * value.commutator_hminus1_norm2;
    value.target_scale = std::pow(graph.z, 1.25L) *
        std::pow(graph.p, 0.75L);
    if (value.target_scale > 0.0L) {
        value.target_ratio = value.upper_envelope /
            value.target_scale;
    }
    value.finite = graph.z > 0.0L && graph.p > 0.0L &&
        std::isfinite(value.upper_envelope) &&
        std::isfinite(value.target_ratio) &&
        graph.hminus1 >= -1e-14L;
    return value;
}

}  // namespace

LocalSldRemainderEnvelopeObjective::
LocalSldRemainderEnvelopeObjective(
    const SpectralDynamics& dynamics,
    TriadSelection selection)
    : dynamics_(dynamics), selection_(selection) {}

LocalSldRemainderEnvelopeValue
LocalSldRemainderEnvelopeObjective::evaluate(
    const SpectralState& state) const {
    return value_from_graph(build_graph(dynamics_, state, selection_));
}

SpectralIncrement
LocalSldRemainderEnvelopeObjective::target_ratio_gradient(
    const SpectralState& state) const {
    const EnvelopeGraph graph = build_graph(
        dynamics_, state, selection_);
    SpectralIncrement result(state.waves.size());
    if (!(graph.z > 0.0L) || !(graph.p > 0.0L)) {
        return result;
    }
    const LocalSldRemainderEnvelopeValue value =
        value_from_graph(graph);

    SpectralIncrement grad_s = graph.ab;
    add_scaled(
        grad_s,
        dynamics_.advection_vjp_direct_partition(
            state, graph.au, selection_),
        1.0L);
    SpectralIncrement grad_j = laplacian_weight(
        state, graph.d, 1);
    add_scaled(
        grad_j,
        d_vjp(
            dynamics_, state, graph, graph.au, selection_),
        1.0L);
    SpectralIncrement grad_h = d_vjp(
        dynamics_, state, graph, graph.inverse_d, selection_);
    for (ComplexVector& mode : grad_h) {
        for (SpectralComplex& component : mode) {
            component *= 2.0L;
        }
    }

    const SpectralReal c_derivative =
        2.0L * graph.c * graph.h3 + graph.j;
    const SpectralReal d_s = graph.s / graph.z +
        3.0L * c_derivative / (4.0L * graph.p);
    const SpectralReal d_z = -graph.s * graph.s /
        (2.0L * graph.z * graph.z);
    const SpectralReal d_p = -c_derivative * graph.c / graph.p;
    add_scaled(result, grad_s, d_s);
    add_scaled(result, graph.au, 2.0L * d_z);
    add_scaled(result, graph.aau, 2.0L * d_p);
    add_scaled(
        result,
        laplacian_weight(state, state.velocity, 3),
        2.0L * graph.c * graph.c);
    add_scaled(result, grad_j, graph.c);
    add_scaled(result, grad_h, 0.25L);

    const SpectralReal factor = 1.0L / value.target_scale;
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= factor;
        }
    }
    add_scaled(result, graph.au,
               -2.5L * value.target_ratio / graph.z);
    add_scaled(result, graph.aau,
               -1.5L * value.target_ratio / graph.p);
    return result;
}

}  // namespace lemma
