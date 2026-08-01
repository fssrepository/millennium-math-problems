#pragma once

#include "local_quartic_closure_objective.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldSignatureBlockReport {
    std::array<SpectralInteger, 3> squared_lengths{1, 1, 2};
    bool equal_low_doubling_family = false;
    std::size_t dominant_interactions = 0;
    std::size_t remainder_interactions = 0;
    SpectralReal advection_reconstruction_error = 0.0L;
    SpectralReal stretching_reconstruction_error = 0.0L;
    SpectralReal bracket_reconstruction_error = 0.0L;
    SpectralReal quotient_reconstruction_error = 0.0L;
    SpectralReal dominant_advection_norm2 = 0.0L;
    SpectralReal remainder_advection_norm2 = 0.0L;
    SpectralReal advection_cross_pairing = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal dominant_stretching = 0.0L;
    SpectralReal remainder_stretching = 0.0L;
    SpectralReal full_bracket = 0.0L;
    SpectralReal dominant_closed_bracket = 0.0L;
    SpectralReal remainder_closed_bracket = 0.0L;
    SpectralReal cross_bracket = 0.0L;
    SpectralReal full_signed_sld_ratio = 0.0L;
    SpectralReal dominant_closed_sld_ratio = 0.0L;
    SpectralReal remainder_closed_sld_ratio = 0.0L;
    SpectralReal cross_sld_ratio = 0.0L;
    SpectralReal dominant_absolute_fraction = 0.0L;
    SpectralReal normalization_initial_frequency = 0.0L;
    SpectralReal normalization_initial_ep_shift = 0.0L;
    int evolved_steps = 0;
    SpectralReal evolved_time = 0.0L;
    SpectralReal viscosity = 0.0L;
    bool frozen_initial_normalization = false;
    LocalQuarticClosureObjectiveValue full;
    LocalQuarticClosureObjectiveValue dominant;
    LocalQuarticClosureObjectiveValue remainder;
    bool exact_decomposition = false;
    bool candidate_block_bound_proved = false;
    bool finite = false;
};

class LocalSldSignatureBlock {
public:
    [[nodiscard]] static LocalSldSignatureBlockReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        std::array<SpectralInteger, 3> squared_lengths = {1, 1, 2},
        bool equal_low_doubling_family = false,
        SpectralReal initial_frequency = 0.0L,
        SpectralReal initial_ep_shift = 0.0L);
};

struct LocalSldSignatureBlockCliOptions {
    std::string state_path;
    std::string certificate_path;
    std::array<SpectralInteger, 3> squared_lengths{1, 1, 2};
    bool equal_low_doubling_family = false;
    int threads = 12;
    int evolve_steps = 0;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
};

class LocalSldSignatureBlockCli {
public:
    [[nodiscard]] static LocalSldSignatureBlockCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldSignatureBlockCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
