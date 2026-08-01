#pragma once

#include "spectral_adjoint.hpp"
#include "initial_sobolev_constraint.hpp"

#include <string>
#include <vector>

namespace lemma {

struct GradientSearchOptions {
    int iterations = 8;
    int line_search_steps = 16;
    int trajectory_steps = 1;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
    SpectralReal initial_step = 0.2L;
    std::string objective = "max-q";
    int sobolev_order = 0;
    SpectralReal sobolev_cap = 0.0L;
};

struct GradientIterationRecord {
    int iteration = 0;
    SpectralReal objective_before = 0.0L;
    SpectralReal objective_after = 0.0L;
    SpectralReal projected_gradient_norm = 0.0L;
    SpectralReal accepted_step = 0.0L;
    SpectralReal sobolev_value = 0.0L;
    int line_search_evaluations = 0;
    bool accepted = false;
};

struct GradientSearchResult {
    SpectralState state;
    SpectralReal initial_objective = 0.0L;
    SpectralReal objective = 0.0L;
    SpectralReal final_projected_gradient_norm = 0.0L;
    SpectralReal final_sobolev_value = 0.0L;
    int objective_step = 0;
    int iterations = 0;
    int accepted_steps = 0;
    int trajectory_evaluations = 0;
    std::vector<GradientIterationRecord> trace;
};

class GradientAdversary {
public:
    GradientAdversary(const SpectralDynamics& dynamics,
                      const SpectralObjective& objective,
                      const SpectralAdjoint& adjoint);

    [[nodiscard]] GradientSearchResult maximize_q(
        const SpectralState& initial,
        const GradientSearchOptions& options) const;

private:
    [[nodiscard]] SpectralReal objective_value(
        const SpectralState& initial,
        const GradientSearchOptions& options) const;

    const SpectralDynamics& dynamics_;
    const SpectralObjective& objective_;
    const SpectralAdjoint& adjoint_;
};

}  // namespace lemma
