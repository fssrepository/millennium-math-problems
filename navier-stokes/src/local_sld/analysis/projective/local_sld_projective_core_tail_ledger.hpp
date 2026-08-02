#pragma once

#include "local_quartic_closure_objective.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

using LocalSldProjectiveCoreSignature =
    std::array<SpectralInteger, 3>;

struct LocalSldProjectiveCoreTailBlock {
    SpectralReal outer_square = 0.0L;
    SpectralReal advected_commutator = 0.0L;
    SpectralReal advecting_nested = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal palinstrophy_normalization = 0.0L;
    SpectralReal bracket = 0.0L;
    SpectralReal power_one = 0.0L;
};

struct LocalSldProjectiveCoreTailReport {
    int cutoff = 0;
    int threads = 1;
    bool excludes_signature_123 = false;
    bool excludes_triple_family = false;
    std::size_t requested_core_shape_count = 0;
    std::size_t active_core_shape_count = 0;
    std::size_t tail_shape_count = 0;
    std::size_t core_interaction_count = 0;
    std::size_t tail_interaction_count = 0;
    std::vector<LocalSldProjectiveCoreSignature> requested_core;
    std::vector<LocalSldProjectiveCoreSignature> active_core;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_local_stretching = 0.0L;
    SpectralReal selected_stretching = 0.0L;
    SpectralReal core_stretching = 0.0L;
    SpectralReal tail_stretching = 0.0L;
    SpectralReal core_palinstrophy_cross = 0.0L;
    SpectralReal tail_palinstrophy_cross = 0.0L;
    SpectralReal power_one_scale = 0.0L;
    SpectralReal selected_bracket = 0.0L;
    SpectralReal selected_power_one = 0.0L;
    LocalSldProjectiveCoreTailBlock core;
    LocalSldProjectiveCoreTailBlock core_tail;
    LocalSldProjectiveCoreTailBlock tail;
    SpectralReal stretching_partition_error = 0.0L;
    SpectralReal palinstrophy_cross_partition_error = 0.0L;
    SpectralReal bracket_partition_error = 0.0L;
    bool exact_core_tail_decomposition = false;
    bool fixed_core_internal_bound_proved = false;
    bool core_tail_bound_proved = false;
    bool growing_tail_internal_bound_proved = false;
    bool full_local_lemma_proved = false;
};

class LocalSldProjectiveCoreTailLedger {
public:
    [[nodiscard]] static LocalSldProjectiveCoreTailReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        std::vector<LocalSldProjectiveCoreSignature> core,
        int threads = 12,
        bool exclude_signature_123 = false,
        bool exclude_triple_family = false);
};

struct LocalSldProjectiveCoreTailCliOptions {
    std::string state_path;
    std::string certificate_path;
    std::vector<LocalSldProjectiveCoreSignature> core;
    SpectralInteger core_maximum_height = 0;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveCoreTailCli {
public:
    [[nodiscard]] static LocalSldProjectiveCoreTailCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveCoreTailCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
