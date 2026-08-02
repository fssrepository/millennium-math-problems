#pragma once

#include "local_sld_response_hierarchy.hpp"

#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace lemma {

struct LocalSldResponseDiagonalOptions {
    std::vector<std::string> state_paths;
    std::string certificate_path;
    int maximum_depth = 16;
    int threads = 12;
    SpectralReal radius = 1.25L;
};

struct LocalSldResponseDiagonalRow {
    std::string state_path;
    int cutoff = 0;
    int safe_order_count = 0;
    int first_truncated_order = -1;
    bool safe_shells_inside_cutoff = true;
    SpectralReal safe_projection_energy = 0.0L;
    SpectralReal safe_weighted_l1 = 0.0L;
    SpectralReal safe_weighted_tail_after_cubic = 0.0L;
    SpectralReal last_safe_coefficient = 0.0L;
    SpectralReal first_truncated_coefficient = 0.0L;
    SpectralReal quadratic_majorant = 0.0L;
    SpectralReal quadratic_majorant_bound = 0.0L;
    SpectralReal quadratic_majorant_ratio = 0.0L;
    SpectralReal common_coefficient_l2_difference = 0.0L;
    SpectralReal common_weighted_l1_difference = 0.0L;
    std::vector<SpectralReal> safe_coefficients;
    std::vector<int> safe_highest_shells;
};

struct LocalSldResponseDiagonalReport {
    int maximum_depth = 0;
    SpectralReal radius = 1.0L;
    SpectralReal maximum_safe_weighted_l1 = 0.0L;
    SpectralReal maximum_common_coefficient_l2_difference = 0.0L;
    SpectralReal maximum_quadratic_majorant_ratio = 0.0L;
    bool finite_cutoff_diagonal_is_not_a_proof = true;
    std::vector<LocalSldResponseDiagonalRow> rows;
};

class LocalSldResponseDiagonal {
public:
    [[nodiscard]] static LocalSldResponseDiagonalReport analyze(
        const SpectralDynamics& dynamics,
        const std::vector<std::pair<std::string, SpectralState>>& states,
        int maximum_depth,
        SpectralReal radius);
};

class LocalSldResponseDiagonalCli {
public:
    [[nodiscard]] static LocalSldResponseDiagonalOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldResponseDiagonalOptions& options,
        std::ostream& out);
};

}  // namespace lemma
