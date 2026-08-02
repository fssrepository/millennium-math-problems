#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationCauchyScanRow {
    SpectralInteger core_maximum_height = 0;
    int cutoff = 0;
    SpectralReal full_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_aggregate_h1_norm2 = 0.0L;
    SpectralReal tail_aggregate_h2_norm2 = 0.0L;
    SpectralReal cauchy_bound_power_one = 0.0L;
    SpectralReal squared_cauchy_bound_power_one = 0.0L;
    SpectralReal quarter_compensated_bound = 0.0L;
    SpectralReal half_compensated_squared_bound = 0.0L;
    SpectralReal numerator = 0.0L;
    SpectralReal denominator = 0.0L;
    std::string state_path;
    bool finite = false;
};

struct LocalSldProjectiveNormalizationCauchyScanReport {
    std::vector<LocalSldProjectiveNormalizationCauchyScanRow> rows;
    SpectralReal fitted_height_slope = 0.0L;
    SpectralReal maximum_bound = 0.0L;
    SpectralReal minimum_bound = 0.0L;
    SpectralReal fitted_quarter_compensated_slope = 0.0L;
    SpectralReal maximum_quarter_compensated_bound = 0.0L;
    bool every_row_finite = false;
    bool finite_scan_is_not_a_proof = true;
    bool uniform_height_decay_proved = false;
};

struct LocalSldProjectiveNormalizationCauchyScanOptions {
    std::vector<SpectralInteger> core_maximum_heights;
    std::vector<std::string> state_paths;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    int threads = 12;
};

class LocalSldProjectiveNormalizationCauchyScanCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveNormalizationCauchyScanOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationCauchyScanOptions& options,
        std::ostream& out);
};

}  // namespace lemma
