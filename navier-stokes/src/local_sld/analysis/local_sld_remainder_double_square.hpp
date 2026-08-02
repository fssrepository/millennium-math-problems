#pragma once

#include "local_quartic_closure_objective.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldRemainderDoubleSquareReport {
    int cutoff = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
    SpectralReal projection_coefficient = 0.0L;
    SpectralReal full_bracket = 0.0L;
    SpectralReal target_scale = 0.0L;
    SpectralReal signed_target_ratio = 0.0L;
    SpectralReal first_square_norm2 = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal projected_h3_correction = 0.0L;
    SpectralReal commutator_pairing = 0.0L;
    SpectralReal projected_commutator_pairing = 0.0L;
    SpectralReal commutator_hminus1_norm2 = 0.0L;
    SpectralReal second_square_norm2 = 0.0L;
    SpectralReal completed_upper_envelope = 0.0L;
    SpectralReal upper_envelope_target_ratio = 0.0L;
    SpectralReal negative_second_square_target_ratio = 0.0L;
    SpectralReal first_completion_error = 0.0L;
    SpectralReal commutator_reconstruction_error = 0.0L;
    SpectralReal stretching_vjp_reconstruction_error = 0.0L;
    SpectralReal second_completion_error = 0.0L;
    bool exact_identity = false;
    bool cutoff_independent_upper_bound_proved = false;
};

class LocalSldRemainderDoubleSquare {
public:
    [[nodiscard]] static LocalSldRemainderDoubleSquareReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state);
};

struct LocalSldRemainderDoubleSquareCliOptions {
    std::string state_path;
    std::string certificate_path;
    int threads = 12;
};

class LocalSldRemainderDoubleSquareCli {
public:
    [[nodiscard]] static LocalSldRemainderDoubleSquareCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldRemainderDoubleSquareCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
