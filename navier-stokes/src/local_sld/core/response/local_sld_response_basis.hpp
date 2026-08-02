#pragma once

#include "spectral_dynamics.hpp"

#include <string>
#include <vector>

namespace lemma {

struct LocalSldResponseBasisElement {
    SpectralState state;
    std::string label;
    int response_order = -1;
    int highest_active_shell = 0;
    int analytic_degree = 0;
    bool scalar_response = false;
};

class LocalSldResponseBasis {
public:
    [[nodiscard]] static std::vector<LocalSldResponseBasisElement> build(
        const SpectralDynamics& dynamics,
        int cutoff,
        int response_depth,
        bool include_transverse_two_one_one = false,
        bool include_three_one_zero_orbits = false);

    [[nodiscard]] static int highest_active_shell(
        const SpectralState& state);
    [[nodiscard]] static SpectralReal maximum_gram_error(
        const std::vector<LocalSldResponseBasisElement>& basis);

    [[nodiscard]] static bool append_orthonormalized(
        const SpectralDynamics& dynamics,
        std::vector<LocalSldResponseBasisElement>& basis,
        SpectralState candidate,
        std::string label,
        int response_order,
        bool scalar_response,
        int analytic_degree = -1);

    [[nodiscard]] static std::vector<LocalSldResponseBasisElement>
    graded_orthonormalize(
        const SpectralDynamics& dynamics,
        std::vector<LocalSldResponseBasisElement> candidates);
};

}  // namespace lemma
