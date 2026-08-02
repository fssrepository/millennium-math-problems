#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationTailRow {
    SpectralInteger core_maximum_height = 0;
    int cutoff = 0;
    SpectralReal open_palinstrophy_normalization_power_one = 0.0L;
    SpectralReal selected_stretching_tail_cross_power_one = 0.0L;
    SpectralReal core_stretching_tail_cross_power_one = 0.0L;
    SpectralReal tail_stretching_core_cross_power_one = 0.0L;
    SpectralReal tail_stretching_tail_cross_power_one = 0.0L;
    SpectralReal joint_cross_tail_absolute_envelope = 0.0L;
    SpectralReal joint_cross_tail_alignment = 0.0L;
    SpectralReal canonical_two_term_absolute_envelope = 0.0L;
    SpectralReal canonical_two_term_alignment = 0.0L;
    SpectralReal selected_stretching_tail_cross_fraction = 0.0L;
    SpectralReal core_stretching_tail_cross_fraction = 0.0L;
    SpectralReal tail_stretching_core_cross_fraction = 0.0L;
    SpectralReal tail_stretching_tail_cross_fraction = 0.0L;
    SpectralReal selected_stretching_tail_cross_cauchy_bound = 0.0L;
    SpectralReal core_stretching_tail_cross_cauchy_bound = 0.0L;
    SpectralReal tail_stretching_core_cross_cauchy_bound = 0.0L;
    SpectralReal tail_stretching_tail_cross_cauchy_bound = 0.0L;
    SpectralReal selected_stretching_tail_cross_cauchy_ratio = 0.0L;
    SpectralReal core_stretching_tail_cross_cauchy_ratio = 0.0L;
    SpectralReal tail_stretching_core_cross_cauchy_ratio = 0.0L;
    SpectralReal tail_stretching_tail_cross_cauchy_ratio = 0.0L;
    SpectralReal joint_cross_tail_cauchy_bound = 0.0L;
    SpectralReal joint_cross_tail_cauchy_ratio = 0.0L;
    SpectralReal core_stretching_h1_alignment = 0.0L;
    SpectralReal tail_stretching_h1_alignment = 0.0L;
    SpectralReal core_palinstrophy_cross_h2_alignment = 0.0L;
    SpectralReal tail_palinstrophy_cross_h2_alignment = 0.0L;
    SpectralReal core_stretching = 0.0L;
    SpectralReal tail_stretching = 0.0L;
    SpectralReal core_palinstrophy_cross = 0.0L;
    SpectralReal tail_palinstrophy_cross = 0.0L;
    SpectralReal factorization_error = 0.0L;
    SpectralReal two_term_factorization_error = 0.0L;
    std::string dominant_channel;
    std::string canonical_dominant_channel;
    std::string state_path;
};

struct LocalSldProjectiveNormalizationTailReport {
    std::vector<LocalSldProjectiveNormalizationTailRow> rows;
    SpectralReal fitted_open_height_slope = 0.0L;
    SpectralReal fitted_selected_stretching_tail_cross_height_slope = 0.0L;
    SpectralReal fitted_core_stretching_tail_cross_height_slope = 0.0L;
    SpectralReal fitted_tail_stretching_core_cross_height_slope = 0.0L;
    SpectralReal fitted_tail_stretching_tail_cross_height_slope = 0.0L;
    SpectralReal fitted_joint_cross_tail_cauchy_bound_height_slope = 0.0L;
    SpectralReal fitted_selected_stretching_tail_cross_cauchy_bound_height_slope = 0.0L;
    SpectralReal maximum_factorization_error = 0.0L;
    SpectralReal minimum_joint_cross_tail_alignment = 0.0L;
    SpectralReal maximum_joint_cross_tail_alignment = 0.0L;
    SpectralReal minimum_joint_cross_tail_cauchy_ratio = 0.0L;
    SpectralReal maximum_joint_cross_tail_cauchy_ratio = 0.0L;
    SpectralReal maximum_individual_cauchy_ratio = 0.0L;
    std::size_t core_stretching_tail_cross_dominance_count = 0;
    std::size_t tail_stretching_core_cross_dominance_count = 0;
    std::size_t tail_stretching_tail_cross_dominance_count = 0;
    std::size_t selected_stretching_tail_cross_canonical_dominance_count = 0;
    std::size_t tail_stretching_core_cross_canonical_dominance_count = 0;
    bool dominant_channel_switch_observed = false;
    bool canonical_dominant_channel_switch_observed = false;
    bool every_matrix_exact = false;
    bool finite_optimized_scan_is_not_a_proof = true;
    bool uniform_joint_cross_tail_bound_proved = false;
};

struct LocalSldProjectiveNormalizationTailCliOptions {
    std::vector<SpectralInteger> core_maximum_heights;
    std::vector<std::string> state_paths;
    std::string certificate_path;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveNormalizationTailScanCli {
public:
    [[nodiscard]] static LocalSldProjectiveNormalizationTailCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationTailCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
