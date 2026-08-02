#pragma once

#include "local_quartic_closure_objective.hpp"
#include "local_sld_projected_square.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldDoublingShellRow {
    int shell = 0;
    SpectralInteger minimum_input_squared = 0;
    SpectralInteger maximum_input_squared_exclusive = 0;
    SpectralReal advection_norm2 = 0.0L;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
};

struct LocalSldDoublingShellMatrixEntry {
    int left_shell = 0;
    int right_shell = 0;
    SpectralReal outer_square = 0.0L;
    SpectralReal advected_commutator = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal palinstrophy_normalization = 0.0L;
    SpectralReal advecting_nested = 0.0L;
    SpectralReal total = 0.0L;
};

struct LocalSldDoublingShellReport {
    int cutoff = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal family_stretching = 0.0L;
    SpectralReal family_palinstrophy_cross = 0.0L;
    SpectralReal absolute_lqc7_ratio = 0.0L;
    SpectralReal signed_local_sld_ratio = 0.0L;
    SpectralReal absolute_target_scale_ratio = 0.0L;
    LocalSldProjectedSquareReport projected_square;
    SpectralReal outer_square_total = 0.0L;
    SpectralReal advected_commutator_total = 0.0L;
    SpectralReal enstrophy_normalization_total = 0.0L;
    SpectralReal palinstrophy_normalization_total = 0.0L;
    SpectralReal advecting_nested_total = 0.0L;
    SpectralReal non_projected_remainder = 0.0L;
    SpectralReal non_projected_target_scale_ratio = 0.0L;
    SpectralReal projected_matrix_reconstruction_error = 0.0L;
    SpectralReal non_projected_reconstruction_error = 0.0L;
    SpectralReal full_family_bracket = 0.0L;
    SpectralReal reconstructed_family_bracket = 0.0L;
    SpectralReal bracket_reconstruction_error = 0.0L;
    SpectralReal advection_reconstruction_error = 0.0L;
    SpectralReal diagonal_total = 0.0L;
    SpectralReal off_diagonal_total = 0.0L;
    SpectralReal absolute_matrix_sum = 0.0L;
    SpectralReal signed_cancellation_fraction = 0.0L;
    SpectralReal maximum_absolute_entry = 0.0L;
    SpectralReal maximum_absolute_far_structural_entry = 0.0L;
    bool exact_reconstruction = false;
    bool structural_entries_are_neighbor_shell_local = false;
    bool signed_cross_shell_bound_proved = false;
    std::vector<LocalSldDoublingShellRow> shells;
    std::vector<LocalSldDoublingShellMatrixEntry> matrix;
};

class LocalSldDoublingShellLedger {
public:
    [[nodiscard]] static LocalSldDoublingShellReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state);
};

struct LocalSldDoublingShellOptions {
    std::string state_path;
    std::string certificate_path;
    int threads = 12;
    int evolve_steps = 0;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
    std::string backend = "auto";
    int two_scale_axis = 0;
    SpectralReal high_to_low_energy_ratio = 0.0L;
    SpectralReal two_scale_response_angle = 0.24L;
    std::string state_output_path;
};

class LocalSldDoublingShellCli {
public:
    [[nodiscard]] static LocalSldDoublingShellOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldDoublingShellOptions& options,
        std::ostream& out);
};

}  // namespace lemma
