#pragma once

#include "helical_sector_adjoint.hpp"

#include <vector>

namespace lemma {

struct HelicalTrajectoryAdversaryOptions {
    HelicalSectorSelection selection =
        HelicalSectorSelection::heterochiral();
    int iterations = 8;
    int line_search_steps = 16;
    int trajectory_steps = 2;
    SpectralReal initial_step = 0.1L;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
};

struct HelicalTrajectoryAdversaryTraceRow {
    int iteration = 0;
    SpectralReal objective = 0.0L;
    SpectralReal projected_gradient_norm = 0.0L;
    SpectralReal accepted_step = 0.0L;
    bool accepted = false;
};

struct HelicalTrajectoryAdversaryResult {
    SpectralState state;
    SpectralReal initial_objective = 0.0L;
    SpectralReal objective = 0.0L;
    int accepted_steps = 0;
    int evaluations = 0;
    std::vector<HelicalTrajectoryAdversaryTraceRow> trace;
};

class HelicalTrajectoryAdversary {
public:
    [[nodiscard]] static HelicalTrajectoryAdversaryResult maximize(
        const SpectralState& initial,
        const HelicalTrajectoryAdversaryOptions& options,
        const HelicalSectorAdjoint& adjoint);
};

}  // namespace lemma
