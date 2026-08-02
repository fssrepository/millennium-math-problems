#pragma once

#include "local_sld_response_hierarchy.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldResponseTensorOptions {
    std::string certificate_path;
    int cutoff = 4;
    int depth = 5;
    int threads = 12;
    SpectralReal input_radius = 1.5L;
    SpectralReal output_radius = 1.25L;
    SpectralReal entry_tolerance = 1e-14L;
    bool include_transverse_two_one_one = false;
    bool include_three_one_zero_orbits = false;
    int closure_extensions = 0;
};

struct LocalSldResponseTensorBasisRow {
    int basis_index = 0;
    int response_order = -1;
    int highest_active_shell = 0;
    int analytic_degree = 0;
    std::string label;
    bool scalar_response = false;
};

struct LocalSldResponseTensorEntry {
    int output_basis_index = 0;
    int output_order = 0;
    int output_analytic_degree = 0;
    std::string output_label;
    SpectralReal coefficient = 0.0L;
};

struct LocalSldResponseTensorShell {
    int shell = 0;
    SpectralReal norm = 0.0L;
    SpectralReal complement_fraction = 0.0L;
    SpectralReal output_weighted_norm = 0.0L;
};

struct LocalSldResponseTensorDegreeOffsetRow {
    int degree_offset = 0;
    std::size_t retained_entries = 0;
    SpectralReal maximum_absolute_coefficient = 0.0L;
    SpectralReal maximum_weighted_entry = 0.0L;
    SpectralReal maximum_pair_block_contribution = 0.0L;
};

struct LocalSldResponseTensorPair {
    int left_basis_index = 0;
    int right_basis_index = 0;
    int left_order = 0;
    int right_order = 0;
    int left_analytic_degree = 0;
    int right_analytic_degree = 0;
    std::string left_label;
    std::string right_label;
    int left_highest_shell = 0;
    int right_highest_shell = 0;
    SpectralReal full_norm = 0.0L;
    SpectralReal projected_norm = 0.0L;
    SpectralReal complement_norm = 0.0L;
    SpectralReal complement_fraction = 0.0L;
    SpectralReal weighted_projected_ratio = 0.0L;
    SpectralReal weighted_projected_degree_block_ratio = 0.0L;
    SpectralReal input_weighted_complement_norm = 0.0L;
    SpectralReal shell_weighted_complement_ratio = 0.0L;
    SpectralReal h1_complement_norm = 0.0L;
    SpectralReal h2_complement_norm = 0.0L;
    SpectralReal input_weighted_h1_complement_norm = 0.0L;
    SpectralReal input_weighted_h2_complement_norm = 0.0L;
    SpectralReal shell_norm_reconstruction_error = 0.0L;
    int maximum_output_degree_excess = 0;
    SpectralReal maximum_degree_excess_coefficient = 0.0L;
    std::vector<LocalSldResponseTensorEntry> entries;
    std::vector<LocalSldResponseTensorShell> complement_shells;
};

struct LocalSldResponseTensorReport {
    int cutoff = 0;
    int depth = 0;
    SpectralReal input_radius = 1.0L;
    SpectralReal output_radius = 1.0L;
    SpectralReal entry_tolerance = 0.0L;
    bool included_transverse_two_one_one = false;
    bool included_three_one_zero_orbits = false;
    std::size_t basis_size = 0;
    int closure_extensions_requested = 0;
    int closure_extensions_constructed = 0;
    SpectralReal maximum_gram_error = 0.0L;
    SpectralReal maximum_projected_bilinear_constant = 0.0L;
    SpectralReal maximum_projected_degree_block_constant = 0.0L;
    SpectralReal axis_pair_candidate_constant = 0.0L;
    SpectralReal projected_constant_over_axis_candidate = 0.0L;
    SpectralReal maximum_complement_norm = 0.0L;
    SpectralReal maximum_complement_fraction = 0.0L;
    SpectralReal maximum_input_weighted_complement_norm = 0.0L;
    SpectralReal maximum_shell_weighted_complement_ratio = 0.0L;
    SpectralReal maximum_input_weighted_h1_complement_norm = 0.0L;
    SpectralReal maximum_input_weighted_h2_complement_norm = 0.0L;
    SpectralReal finite_combined_weighted_bound = 0.0L;
    SpectralReal maximum_norm_reconstruction_error = 0.0L;
    SpectralReal maximum_shell_norm_reconstruction_error = 0.0L;
    int maximum_output_degree_excess = 0;
    SpectralReal maximum_degree_excess_coefficient = 0.0L;
    std::size_t retained_tensor_entries = 0;
    bool boundary_free_depth = false;
    bool axis_pair_candidate_survives = false;
    bool graded_support_closed = true;
    bool finite_tensor_is_not_a_proof = true;
    std::vector<LocalSldResponseTensorBasisRow> basis;
    std::vector<LocalSldResponseTensorDegreeOffsetRow>
        degree_offset_envelope;
    std::vector<LocalSldResponseTensorPair> pairs;
};

class LocalSldResponseTensor {
public:
    [[nodiscard]] static LocalSldResponseTensorReport analyze(
        const SpectralDynamics& dynamics,
        int cutoff,
        int depth,
        SpectralReal input_radius,
        SpectralReal output_radius,
        SpectralReal entry_tolerance,
        bool include_transverse_two_one_one = false,
        bool include_three_one_zero_orbits = false,
        int closure_extensions = 0);
};

class LocalSldResponseTensorCli {
public:
    [[nodiscard]] static LocalSldResponseTensorOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldResponseTensorOptions& options,
        std::ostream& out);
};

}  // namespace lemma
