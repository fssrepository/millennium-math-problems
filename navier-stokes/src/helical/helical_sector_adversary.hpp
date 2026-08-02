#pragma once

#include "helical_sector_objective.hpp"

#include <vector>

namespace lemma {

struct HelicalSectorAdversaryOptions {
    HelicalSectorSelection selection =
        HelicalSectorSelection::heterochiral();
    int iterations = 8;
    int line_search_steps = 16;
    SpectralReal initial_step = 0.1L;
};

struct HelicalSectorAdversaryTraceRow {
    int iteration = 0;
    SpectralReal objective = 0.0L;
    SpectralReal projected_gradient_norm = 0.0L;
    SpectralReal accepted_step = 0.0L;
    bool accepted = false;
};

struct HelicalSectorAdversaryResult {
    SpectralState state;
    SpectralReal initial_objective = 0.0L;
    SpectralReal objective = 0.0L;
    int accepted_steps = 0;
    int evaluations = 0;
    std::vector<HelicalSectorAdversaryTraceRow> trace;
};

class HelicalSectorAdversary {
public:
    [[nodiscard]] static HelicalSectorAdversaryResult maximize(
        const SpectralState& initial,
        const HelicalSectorAdversaryOptions& options);
};

}  // namespace lemma
