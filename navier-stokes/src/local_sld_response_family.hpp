#pragma once

#include "local_sld_response_hierarchy.hpp"

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace lemma {

struct LocalSldResponseFamilyOptions {
    std::vector<std::string> state_paths;
    std::string certificate_path;
    int depth = 16;
    int threads = 12;
    bool include_transverse_two_one_one = false;
    bool include_three_one_zero_orbits = false;
};

struct LocalSldResponseFamilyRow {
    std::string state_path;
    int cutoff = 0;
    SpectralReal axis_response_energy = 0.0L;
    SpectralReal cubic_response_energy = 0.0L;
    SpectralReal finite_response_energy = 0.0L;
    SpectralReal orbit_energy_increment = 0.0L;
    SpectralReal final_projection_energy = 0.0L;
    SpectralReal final_projection_residual = 0.0L;
    SpectralReal coefficient_l2_difference = 0.0L;
    std::vector<std::string> labels;
    std::vector<SpectralReal> coefficients;
    std::vector<SpectralReal> coefficient_energies;
};

struct LocalSldResponseFamilyReport {
    int depth = 0;
    SpectralReal maximum_coefficient_l2_difference = 0.0L;
    SpectralReal maximum_projection_residual = 0.0L;
    std::vector<LocalSldResponseFamilyRow> rows;
    bool finite_cutoff_family_is_not_a_proof = true;
};

class LocalSldResponseFamily {
public:
    [[nodiscard]] static LocalSldResponseFamilyReport analyze(
        const SpectralDynamics& dynamics,
        const std::vector<std::pair<std::string, SpectralState>>& states,
        int depth,
        bool include_transverse_two_one_one,
        bool include_three_one_zero_orbits);
};

class LocalSldResponseFamilyCli {
public:
    [[nodiscard]] static LocalSldResponseFamilyOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldResponseFamilyOptions& options,
        std::ostream& out);
};

}  // namespace lemma
