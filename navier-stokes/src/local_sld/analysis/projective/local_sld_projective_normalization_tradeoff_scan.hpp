#pragma once

#include "spectral_dynamics.hpp"
#include "spectral_state.hpp"
#include "triad_partition.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationTradeoffRow {
    SpectralReal parameter = 0.0L;
    SpectralReal open_palinstrophy_normalization_power_one = 0.0L;
    SpectralReal squared_open_palinstrophy_normalization_power_one = 0.0L;
    SpectralReal tail_stretching_alignment_squared = 0.0L;
    SpectralReal selected_stretching_h1_alignment_squared = 0.0L;
    SpectralReal tail_palinstrophy_cross_h2_alignment_squared = 0.0L;
    SpectralReal normalization_alignment_product_squared = 0.0L;
    SpectralReal selected_stretching_tail_cross_power_one = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal palinstrophy_over_enstrophy_squared = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal normalization_scale = 0.0L;
    SpectralReal unscaled_open_normalization = 0.0L;
};

struct LocalSldProjectiveNormalizationTradeoffReport {
    int cutoff = 0;
    SpectralInteger core_maximum_height = 0;
    SpectralReal maximum_open_power_one = 0.0L;
    SpectralReal maximum_open_parameter = 0.0L;
    SpectralReal maximum_tail_alignment_squared = 0.0L;
    SpectralReal maximum_tail_alignment_parameter = 0.0L;
    SpectralReal maximum_normalization_alignment_product_squared = 0.0L;
    SpectralReal maximum_normalization_alignment_parameter = 0.0L;
    SpectralReal maximum_selected_channel_power_one = 0.0L;
    SpectralReal maximum_selected_channel_parameter = 0.0L;
    SpectralReal maximum_roughness = 0.0L;
    SpectralReal maximum_roughness_parameter = 0.0L;
    bool every_sample_finite = false;
    bool finite_line_scan_is_not_a_proof = true;
    bool uniform_scale_alignment_tradeoff_proved = false;
    std::vector<LocalSldProjectiveNormalizationTradeoffRow> rows;
};

class LocalSldProjectiveNormalizationTradeoffScan {
public:
    [[nodiscard]] static LocalSldProjectiveNormalizationTradeoffReport
    analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& normalization_state,
        const SpectralState& alignment_state,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        SpectralReal minimum_parameter,
        SpectralReal maximum_parameter,
        int samples,
        int threads);
};

struct LocalSldProjectiveNormalizationTradeoffScanOptions {
    std::string normalization_state_path;
    std::string alignment_state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    SpectralInteger core_maximum_height = 8;
    SpectralReal minimum_parameter = 0.0L;
    SpectralReal maximum_parameter = 1.0L;
    int samples = 25;
    int threads = 12;
};

class LocalSldProjectiveNormalizationTradeoffScanCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveNormalizationTradeoffScanOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationTradeoffScanOptions& options,
        std::ostream& out);
};

}  // namespace lemma
