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
};

struct LocalSldResponseTensorEntry {
    int output_order = 0;
    SpectralReal coefficient = 0.0L;
};

struct LocalSldResponseTensorPair {
    int left_order = 0;
    int right_order = 0;
    int left_highest_shell = 0;
    int right_highest_shell = 0;
    SpectralReal full_norm = 0.0L;
    SpectralReal projected_norm = 0.0L;
    SpectralReal complement_norm = 0.0L;
    SpectralReal complement_fraction = 0.0L;
    SpectralReal weighted_projected_ratio = 0.0L;
    SpectralReal input_weighted_complement_norm = 0.0L;
    std::vector<LocalSldResponseTensorEntry> entries;
};

struct LocalSldResponseTensorReport {
    int cutoff = 0;
    int depth = 0;
    SpectralReal input_radius = 1.0L;
    SpectralReal output_radius = 1.0L;
    SpectralReal entry_tolerance = 0.0L;
    SpectralReal maximum_gram_error = 0.0L;
    SpectralReal maximum_projected_bilinear_constant = 0.0L;
    SpectralReal maximum_complement_norm = 0.0L;
    SpectralReal maximum_complement_fraction = 0.0L;
    SpectralReal maximum_input_weighted_complement_norm = 0.0L;
    SpectralReal maximum_norm_reconstruction_error = 0.0L;
    std::size_t retained_tensor_entries = 0;
    bool boundary_free_depth = false;
    bool finite_tensor_is_not_a_proof = true;
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
        SpectralReal entry_tolerance);
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
