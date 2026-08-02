#pragma once

#include "spectral_dynamics.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldResponseHierarchyOptions {
    std::string state_path;
    std::string certificate_path;
    std::string residual_state_path;
    std::string projected_state_path;
    int depth = 6;
    int threads = 12;
    bool include_transverse_two_one_one = false;
    bool include_three_one_zero_orbits = false;
};

struct LocalSldResponseHierarchyRow {
    int order = 0;
    std::string label;
    SpectralReal coefficient = 0.0L;
    SpectralReal coefficient_energy = 0.0L;
    SpectralReal cumulative_projection_energy = 0.0L;
    SpectralReal projection_residual = 0.0L;
    int highest_active_shell = 0;
};

struct LocalSldResponseHierarchyReport {
    SpectralState residual_state;
    SpectralState projected_state;
    int cutoff = 0;
    int requested_depth = 0;
    int constructed_depth = 0;
    bool included_transverse_two_one_one = false;
    bool included_three_one_zero_orbits = false;
    SpectralReal reference_energy = 0.0L;
    SpectralReal maximum_gram_error = 0.0L;
    SpectralReal final_projection_energy = 0.0L;
    SpectralReal final_projection_residual = 0.0L;
    std::vector<LocalSldResponseHierarchyRow> rows;
    bool finite_search_is_not_a_proof = true;
};

class LocalSldResponseHierarchy {
public:
    [[nodiscard]] static std::vector<SpectralState> build(
        const SpectralDynamics& dynamics,
        int cutoff,
        int depth);
    [[nodiscard]] static LocalSldResponseHierarchyReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& reference,
        int depth,
        bool include_transverse_two_one_one = false,
        bool include_three_one_zero_orbits = false);
};

class LocalSldResponseHierarchyCli {
public:
    [[nodiscard]] static LocalSldResponseHierarchyOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldResponseHierarchyOptions& options,
        std::ostream& out);
};

}  // namespace lemma
