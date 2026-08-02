#include "local_sld_response_basis.hpp"

#include "local_sld_cyclic_basis.hpp"
#include "local_sld_cyclic_orbit_basis.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace lemma {
namespace {

void add_increment(
    SpectralIncrement& target,
    const SpectralIncrement& source) {
    if (target.size() != source.size()) {
        throw std::invalid_argument(
            "response-basis increment layout mismatch");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += source[mode][component];
        }
    }
}

void subtract_projection(
    SpectralState& candidate,
    const SpectralState& previous) {
    const SpectralReal coefficient = LocalSldCyclicBasis::pairing(
        previous.velocity, candidate.velocity);
    for (std::size_t mode = 0; mode < candidate.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            candidate.velocity[mode][component] -=
                coefficient * previous.velocity[mode][component];
        }
    }
}

}  // namespace

std::vector<LocalSldResponseBasisElement> LocalSldResponseBasis::build(
    const SpectralDynamics& dynamics,
    int cutoff,
    int response_depth,
    bool include_transverse_two_one_one,
    bool include_three_one_zero_orbits) {
    if (cutoff < 2 || cutoff > 12 || response_depth < 2 ||
        response_depth > 16) {
        throw std::invalid_argument(
            "response basis requires cutoff 2..12 and depth 2..16");
    }
    if (include_three_one_zero_orbits && cutoff < 3) {
        throw std::invalid_argument(
            "the (3,1,0) response orbits require cutoff at least 3");
    }
    std::vector<LocalSldResponseBasisElement> basis;
    basis.reserve(static_cast<std::size_t>(
        response_depth + (include_transverse_two_one_one ? 1 : 0) +
        (include_three_one_zero_orbits ? 2 : 0)));
    if (!append_orthonormalized(
        dynamics, basis,
        LocalSldCyclicBasis::axis_state(cutoff),
        "response-order-0", 0, true)) {
        throw std::runtime_error("cyclic axis basis state vanished");
    }
    if (!append_orthonormalized(
        dynamics, basis,
        LocalSldCyclicBasis::response_state(
            dynamics, basis.front().state),
        "response-order-1", 1, true)) {
        throw std::runtime_error("first cyclic response state vanished");
    }
    for (int order = 2; order < response_depth; ++order) {
        SpectralState candidate = basis.front().state;
        for (ComplexVector& value : candidate.velocity) {
            value = {};
        }
        for (int left = 0; left < order; ++left) {
            const int right = order - 1 - left;
            add_increment(
                candidate.velocity,
                dynamics.advection_bilinear_direct_partition(
                    candidate,
                    basis[static_cast<std::size_t>(left)].state.velocity,
                    basis[static_cast<std::size_t>(right)].state.velocity,
                    TriadPartition::all));
        }
        if (!append_orthonormalized(
            dynamics, basis, std::move(candidate),
            "response-order-" + std::to_string(order),
            order, true)) {
            throw std::runtime_error(
                "response recursion vanished at order " +
                std::to_string(order));
        }
    }
    if (!include_transverse_two_one_one &&
        !include_three_one_zero_orbits) {
        return basis;
    }
    std::vector<LocalSldResponseBasisElement> candidates = basis;
    auto append_orbit_candidate = [&](SpectralState state,
                                      std::string label) {
        const int highest_shell = highest_active_shell(state);
        candidates.push_back(LocalSldResponseBasisElement{
            std::move(state), std::move(label), -1, highest_shell,
            std::max(0, highest_shell - 1), false});
    };
    if (include_transverse_two_one_one) {
        append_orbit_candidate(
            LocalSldCyclicOrbitBasis::transverse_two_one_one(cutoff),
            "transverse-(2,1,1)-orbit");
    }
    if (include_three_one_zero_orbits) {
        append_orbit_candidate(
            LocalSldCyclicOrbitBasis::forward_three_one_zero(cutoff),
            "forward-(3,1,0)-orbit");
        append_orbit_candidate(
            LocalSldCyclicOrbitBasis::backward_three_one_zero(cutoff),
            "backward-(3,1,0)-orbit");
    }
    return graded_orthonormalize(dynamics, std::move(candidates));
}

int LocalSldResponseBasis::highest_active_shell(
    const SpectralState& state) {
    SpectralReal maximum_mode_energy = 0.0L;
    for (const ComplexVector& value : state.velocity) {
        maximum_mode_energy = std::max(
            maximum_mode_energy,
            std::real(dot_hermitian(value, value)));
    }
    const SpectralReal threshold = 1e-24L * maximum_mode_energy;
    int shell = 0;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        if (std::real(dot_hermitian(
                state.velocity[mode], state.velocity[mode])) <= threshold) {
            continue;
        }
        const WaveVector wave = state.waves[mode];
        shell = std::max(shell, std::max(
            {std::abs(wave.x), std::abs(wave.y), std::abs(wave.z)}));
    }
    return shell;
}

SpectralReal LocalSldResponseBasis::maximum_gram_error(
    const std::vector<LocalSldResponseBasisElement>& basis) {
    SpectralReal error = 0.0L;
    for (std::size_t left = 0; left < basis.size(); ++left) {
        for (std::size_t right = 0; right < basis.size(); ++right) {
            const SpectralReal expected = left == right ? 1.0L : 0.0L;
            error = std::max(error, std::abs(
                LocalSldCyclicBasis::pairing(
                    basis[left].state.velocity,
                    basis[right].state.velocity) - expected));
        }
    }
    return error;
}

bool LocalSldResponseBasis::append_orthonormalized(
    const SpectralDynamics& dynamics,
    std::vector<LocalSldResponseBasisElement>& basis,
    SpectralState candidate,
    std::string label,
    int response_order,
    bool scalar_response,
    int analytic_degree) {
    for (const LocalSldResponseBasisElement& previous : basis) {
        subtract_projection(candidate, previous.state);
    }
    dynamics.enforce_constraints(candidate);
    if (!(SpectralStateOps::energy(candidate) > 1e-28L)) {
        return false;
    }
    SpectralStateOps::normalize_energy(candidate);
    const int highest_shell = highest_active_shell(candidate);
    const int resolved_degree = analytic_degree >= 0
        ? analytic_degree
        : (scalar_response
            ? response_order
            : std::max(0, highest_shell - 1));
    basis.push_back(LocalSldResponseBasisElement{
        std::move(candidate),
        std::move(label),
        response_order,
        highest_shell,
        resolved_degree,
        scalar_response});
    return true;
}

std::vector<LocalSldResponseBasisElement>
LocalSldResponseBasis::graded_orthonormalize(
    const SpectralDynamics& dynamics,
    std::vector<LocalSldResponseBasisElement> candidates) {
    std::stable_sort(
        candidates.begin(), candidates.end(),
        [](const LocalSldResponseBasisElement& left,
           const LocalSldResponseBasisElement& right) {
            return left.analytic_degree < right.analytic_degree;
        });
    std::vector<LocalSldResponseBasisElement> graded_basis;
    graded_basis.reserve(candidates.size());
    for (LocalSldResponseBasisElement& candidate : candidates) {
        const std::string label = candidate.label;
        if (!append_orthonormalized(
                dynamics, graded_basis, std::move(candidate.state),
                label, candidate.response_order,
                candidate.scalar_response,
                candidate.analytic_degree)) {
            throw std::runtime_error(
                "graded response-basis candidate vanished: " + label);
        }
    }
    return graded_basis;
}

}  // namespace lemma
