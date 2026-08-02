#pragma once

#include "local_quartic_closure_objective.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightMatrixEntry {
    int first_shell = 0;
    int second_shell = 0;
    SpectralInteger first_minimum_height = 0;
    SpectralInteger first_maximum_height = 0;
    SpectralInteger second_minimum_height = 0;
    SpectralInteger second_maximum_height = 0;
    SpectralReal outer_square = 0.0L;
    SpectralReal advected_commutator = 0.0L;
    SpectralReal advecting_nested = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal palinstrophy_normalization = 0.0L;
    SpectralReal bracket = 0.0L;
    SpectralReal power_one = 0.0L;
    SpectralReal absolute_component_power_one_envelope = 0.0L;
    SpectralReal commutator_paired_power_one_envelope = 0.0L;
    SpectralReal absolute_fraction = 0.0L;
};

struct LocalSldProjectiveHeightShellSummary {
    int shell = 0;
    SpectralInteger minimum_height = 0;
    SpectralInteger maximum_height = 0;
    std::size_t shape_count = 0;
    std::size_t interaction_count = 0;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
    SpectralReal aggregate_l2_norm2 = 0.0L;
    SpectralReal aggregate_h1_norm2 = 0.0L;
    SpectralReal aggregate_h2_norm2 = 0.0L;
    SpectralReal square_function_l2_norm2 = 0.0L;
    SpectralReal square_function_h1_norm2 = 0.0L;
    SpectralReal square_function_h2_norm2 = 0.0L;
    SpectralReal l2_synthesis_ratio = 0.0L;
    SpectralReal h1_synthesis_ratio = 0.0L;
    SpectralReal h2_synthesis_ratio = 0.0L;
    SpectralReal stretching_alignment_squared = 0.0L;
    SpectralReal stretching_h1_alignment_squared = 0.0L;
    SpectralReal h1_synthesis_stretching_product = 0.0L;
    SpectralReal palinstrophy_cross_alignment_squared = 0.0L;
};

struct LocalSldProjectiveHeightMatrixReport {
    int cutoff = 0;
    int threads = 1;
    bool excludes_signature_123 = false;
    bool excludes_triple_family = false;
    SpectralReal selected_bracket = 0.0L;
    SpectralReal selected_power_one = 0.0L;
    SpectralReal power_one_scale = 0.0L;
    SpectralReal reconstructed_bracket = 0.0L;
    SpectralReal reconstructed_power_one = 0.0L;
    SpectralReal absolute_power_one_sum = 0.0L;
    SpectralReal effective_height_pairs = 0.0L;
    SpectralReal dominant_height_pair_fraction = 0.0L;
    SpectralReal signed_height_pair_alignment = 0.0L;
    SpectralReal bracket_reconstruction_error = 0.0L;
    bool exact_height_matrix_decomposition = false;
    bool finite_height_matrix_is_not_a_proof = true;
    std::vector<LocalSldProjectiveHeightShellSummary> shells;
    std::vector<LocalSldProjectiveHeightMatrixEntry> entries;
};

class LocalSldProjectiveHeightMatrix {
public:
    [[nodiscard]] static LocalSldProjectiveHeightMatrixReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        int threads = 12,
        bool exclude_signature_123 = false,
        bool exclude_triple_family = false);
};

struct LocalSldProjectiveHeightMatrixCliOptions {
    std::string state_path;
    std::string certificate_path;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveHeightMatrixCli {
public:
    [[nodiscard]] static LocalSldProjectiveHeightMatrixCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightMatrixCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
