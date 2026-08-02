#pragma once

#include "local_sld_projective_coherence_ledger.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveFanScanRow {
    int cutoff = 0;
    std::size_t coherent_pairs = 0;
    std::size_t active_positive_modes = 0;
    std::size_t projective_shapes = 0;
    SpectralReal synthesis_ratio = 0.0L;
    SpectralReal synthesis_amplification = 0.0L;
    SpectralReal maximum_output_synthesis_ratio = 0.0L;
    WaveVector maximum_output_wave{};
    SpectralReal bracket_constant_ratio = 0.0L;
    SpectralReal normalized_full_stretching = 0.0L;
    SpectralReal power_one_product = 0.0L;
    SpectralReal projective_reconstruction_error = 0.0L;
    std::string state_path;
};

struct LocalSldProjectiveFanScanReport {
    int threads = 1;
    SpectralReal synthesis_ratio_cutoff_slope = 0.0L;
    SpectralReal power_one_cutoff_slope = 0.0L;
    bool every_projective_reconstruction_exact = false;
    bool finite_scan_is_not_a_proof = true;
    std::vector<LocalSldProjectiveFanScanRow> rows;
};

struct LocalSldProjectiveFanScanOptions {
    int minimum_cutoff = 3;
    int maximum_cutoff = 12;
    int threads = 12;
    std::string certificate_path;
    std::string state_directory;
};

class LocalSldProjectiveFanScan {
public:
    [[nodiscard]] static LocalSldProjectiveFanScanReport analyze(
        const LocalSldProjectiveFanScanOptions& options);
    [[nodiscard]] static LocalSldProjectiveFanScanOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveFanScanOptions& options,
        std::ostream& out);
};

}  // namespace lemma
