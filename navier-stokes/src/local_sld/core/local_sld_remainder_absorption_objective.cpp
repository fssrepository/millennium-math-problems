#include "local_sld_remainder_absorption_objective.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lemma {
namespace {

SpectralReal pairing(const SpectralIncrement& left,
                     const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "remainder absorption layout mismatch");
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
    int power = 1) {
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        SpectralReal weight = 1.0L;
        for (int index = 0; index < power; ++index) {
            weight *= wave2;
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
            "remainder absorption gradient mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += scale * source[mode][component];
        }
    }
}

struct SquareGraph {
    SpectralIncrement au;
    SpectralIncrement aau;
    SpectralIncrement b;
    SpectralIncrement ab;
    SpectralIncrement w;
    SpectralIncrement aw;
    SpectralReal z = 0.0L;
    SpectralReal p = 0.0L;
    SpectralReal s = 0.0L;
    SpectralReal c = 0.0L;
    SpectralReal square = 0.0L;
};

SquareGraph build_square_graph(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    TriadSelection selection) {
    SquareGraph graph;
    graph.au = laplacian_weight(state, state.velocity);
    graph.aau = laplacian_weight(state, graph.au);
    graph.b = dynamics.advection_direct_partition(state, selection);
    graph.ab = laplacian_weight(state, graph.b);
    graph.z = pairing(state.velocity, graph.au);
    graph.p = pairing(graph.au, graph.au);
    graph.s = pairing(graph.au, graph.b);
    graph.c = 3.0L * graph.s / (4.0L * graph.p);
    graph.w = graph.b;
    add_scaled(graph.w, graph.au, -graph.c);
    graph.aw = laplacian_weight(state, graph.w);
    graph.square = pairing(graph.w, graph.aw);
    return graph;
}

SpectralIncrement first_square_gradient(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const SquareGraph& graph,
    TriadSelection selection) {
    SpectralIncrement result =
        dynamics.advection_vjp_direct_partition(
            state, graph.aw, selection);
    for (ComplexVector& mode : result) {
        for (SpectralComplex& component : mode) {
            component *= 2.0L;
        }
    }
    add_scaled(
        result,
        laplacian_weight(state, graph.aw),
        -2.0L * graph.c);

    SpectralIncrement grad_s = graph.ab;
    add_scaled(
        grad_s,
        dynamics.advection_vjp_direct_partition(
            state, graph.au, selection),
        1.0L);
    const SpectralReal dc_cotangent =
        -2.0L * pairing(graph.aw, graph.au);
    add_scaled(
        result, grad_s,
        dc_cotangent * 3.0L / (4.0L * graph.p));
    add_scaled(
        result, graph.aau,
        dc_cotangent * (-2.0L * graph.c / graph.p));
    return result;
}

}  // namespace

LocalSldRemainderAbsorptionObjective::
LocalSldRemainderAbsorptionObjective(
    const SpectralDynamics& dynamics,
    SpectralReal theta,
    TriadSelection selection)
    : dynamics_(dynamics), theta_(theta), selection_(selection) {
    if (!(theta_ >= 0.0L) || !(theta_ <= 1.0L) ||
        !std::isfinite(theta_)) {
        throw std::invalid_argument(
            "remainder absorption theta must be in [0,1]");
    }
}

LocalSldRemainderAbsorptionValue
LocalSldRemainderAbsorptionObjective::evaluate(
    const SpectralState& state) const {
    const LocalQuarticClosureObjectiveValue closure =
        LocalQuarticClosureObjective(
            dynamics_, selection_).evaluate(state);
    const SquareGraph graph = build_square_graph(
        dynamics_, state, selection_);
    LocalSldRemainderAbsorptionValue value;
    value.theta = theta_;
    value.signed_bracket = closure.signed_two_entry_bracket;
    value.first_square_norm2 = graph.square;
    value.target_scale = closure.lqc3_target_scale;
    value.signed_bracket_ratio = closure.signed_lqc3_target_ratio;
    if (value.target_scale > 0.0L) {
        value.first_square_ratio = graph.square / value.target_scale;
        value.absorption_ratio = value.signed_bracket_ratio +
            (1.0L - theta_) * value.first_square_ratio;
    }
    value.finite = closure.finite && graph.square >= -1e-14L &&
        std::isfinite(value.absorption_ratio);
    return value;
}

SpectralIncrement
LocalSldRemainderAbsorptionObjective::absorption_ratio_gradient(
    const SpectralState& state) const {
    const LocalQuarticClosureObjective closure(
        dynamics_, selection_);
    SpectralIncrement result =
        closure.signed_lqc3_target_ratio_gradient(state);
    if (theta_ == 1.0L) {
        return result;
    }
    const LocalQuarticClosureObjectiveValue closure_value =
        closure.evaluate(state);
    const SquareGraph graph = build_square_graph(
        dynamics_, state, selection_);
    SpectralIncrement square_gradient = first_square_gradient(
        dynamics_, state, graph, selection_);
    const SpectralReal factor = 1.0L /
        closure_value.lqc3_target_scale;
    for (ComplexVector& mode : square_gradient) {
        for (SpectralComplex& component : mode) {
            component *= factor;
        }
    }
    const SpectralReal square_ratio = graph.square * factor;
    add_scaled(
        square_gradient, graph.au,
        -2.5L * square_ratio / graph.z);
    add_scaled(
        square_gradient, graph.aau,
        -1.5L * square_ratio / graph.p);
    add_scaled(result, square_gradient, 1.0L - theta_);
    return result;
}

}  // namespace lemma
