#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationSatelliteRow {
    int cutoff = 0;
    SpectralReal satellite_amplitude = 0.0L;
    SpectralReal satellite_energy_fraction = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal selected_aggregate_h1_norm2 = 0.0L;
    SpectralReal tail_aggregate_h2_norm2 = 0.0L;
    SpectralReal cauchy_bound_power_one = 0.0L;
    SpectralReal selected_channel_power_one = 0.0L;
    SpectralReal actual_over_cauchy_bound = 0.0L;
    SpectralReal quarter_compensated_bound = 0.0L;
    std::string state_path;
    bool finite = false;
};

struct LocalSldProjectiveNormalizationSatelliteReport {
    std::vector<LocalSldProjectiveNormalizationSatelliteRow> rows;
    SpectralReal fitted_cauchy_cutoff_slope = 0.0L;
    SpectralReal fitted_actual_cutoff_slope = 0.0L;
    SpectralReal maximum_cauchy_bound = 0.0L;
    SpectralReal maximum_selected_channel = 0.0L;
    bool every_row_finite = false;
    bool finite_satellite_scan_is_not_a_proof = true;
    bool uniform_cauchy_majorant_proved = false;
};

struct LocalSldProjectiveNormalizationSatelliteOptions {
    std::string base_state_path;
    std::string state_directory;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    std::vector<int> cutoffs;
    int base_cutoff = 3;
    SpectralInteger core_maximum_height = 8;
    SpectralReal satellite_coefficient = 1.0L;
    int satellite_mode_pairs = 1;
    std::uint64_t seed = 20260802;
    int threads = 12;
};

class LocalSldProjectiveNormalizationSatelliteCli {
public:
    [[nodiscard]] static LocalSldProjectiveNormalizationSatelliteOptions
    parse(int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationSatelliteOptions& options,
        std::ostream& out);
};

}  // namespace lemma
