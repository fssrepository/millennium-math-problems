#pragma once

#include "local_sld_projected_square.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldDoublingScaleRow {
    int scale = 0;
    SpectralReal response_angle = 0.0L;
    SpectralReal high_to_low_energy_ratio = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal bracket = 0.0L;
    SpectralReal signed_target_ratio = 0.0L;
    SpectralReal full_local_bracket = 0.0L;
    SpectralReal full_local_target_ratio = 0.0L;
    SpectralReal remainder_closed_bracket = 0.0L;
    SpectralReal remainder_closed_target_ratio = 0.0L;
    SpectralReal mixed_bracket = 0.0L;
    SpectralReal mixed_target_ratio = 0.0L;
    SpectralReal projected = 0.0L;
    SpectralReal projected_target_ratio = 0.0L;
    SpectralReal non_projected = 0.0L;
    SpectralReal non_projected_target_ratio = 0.0L;
    SpectralReal cancellation_fraction = 0.0L;
    SpectralReal square_identity_error = 0.0L;
};

struct LocalSldDoublingScaleScanReport {
    int minimum_scale = 0;
    int maximum_scale = 0;
    int angle_count = 0;
    SpectralReal energy_decay_power = 0.0L;
    LocalSldDoublingScaleRow maximum_signed_row;
    LocalSldDoublingScaleRow maximum_absolute_row;
    LocalSldDoublingScaleRow maximum_non_projected_row;
    LocalSldDoublingScaleRow maximum_remainder_row;
    LocalSldDoublingScaleRow maximum_mixed_row;
    std::vector<LocalSldDoublingScaleRow> rows;
    bool every_square_identity_verified = false;
    bool every_block_decomposition_verified = false;
    bool cutoff_independent_bound_proved = false;
};

struct LocalSldDoublingScaleScanOptions {
    int minimum_scale = 2;
    int maximum_scale = 12;
    SpectralReal minimum_angle = -1.2L;
    SpectralReal maximum_angle = 1.2L;
    int angle_count = 25;
    SpectralReal energy_decay_power = 2.75L;
    int threads = 12;
    std::string certificate_path;
};

class LocalSldDoublingScaleScan {
public:
    [[nodiscard]] static LocalSldDoublingScaleScanReport analyze(
        const LocalSldDoublingScaleScanOptions& options);
    [[nodiscard]] static LocalSldDoublingScaleScanOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldDoublingScaleScanOptions& options,
        std::ostream& out);
};

}  // namespace lemma
