#pragma once

#include "local_sld_projective_core_tail_ledger.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveCoreTailScanRow {
    int cutoff = 0;
    std::size_t active_core_shape_count = 0;
    std::size_t tail_shape_count = 0;
    SpectralReal selected_power_one = 0.0L;
    SpectralReal core_power_one = 0.0L;
    SpectralReal core_tail_power_one = 0.0L;
    SpectralReal tail_power_one = 0.0L;
    SpectralReal open_signed_power_one = 0.0L;
    SpectralReal open_absolute_power_one = 0.0L;
    SpectralReal open_signed_fraction = 0.0L;
    SpectralReal open_absolute_fraction = 0.0L;
    SpectralReal decomposition_error = 0.0L;
};

struct LocalSldProjectiveCoreTailScanReport {
    std::vector<LocalSldProjectiveCoreTailScanRow> rows;
    SpectralReal fitted_open_signed_cutoff_slope = 0.0L;
    SpectralReal fitted_open_absolute_cutoff_slope = 0.0L;
    int full_active_core_start_cutoff = 0;
    SpectralReal fitted_full_active_core_open_signed_cutoff_slope = 0.0L;
    SpectralReal fitted_full_active_core_open_absolute_cutoff_slope = 0.0L;
    bool every_decomposition_exact = false;
    bool fixed_core_internal_bound_proved = false;
    bool uniform_open_bound_proved = false;
    bool full_local_lemma_proved = false;
};

class LocalSldProjectiveCoreTailScan {
public:
    [[nodiscard]] static LocalSldProjectiveCoreTailScanReport analyze(
        const SpectralDynamics& dynamics,
        const std::vector<SpectralState>& states,
        const std::vector<LocalSldProjectiveCoreSignature>& core,
        int threads = 12,
        bool exclude_signature_123 = false,
        bool exclude_triple_family = false);
};

struct LocalSldProjectiveCoreTailScanCliOptions {
    std::vector<std::string> state_paths;
    std::vector<LocalSldProjectiveCoreSignature> core;
    SpectralInteger core_maximum_height = 0;
    std::string certificate_path;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveCoreTailScanCli {
public:
    [[nodiscard]] static LocalSldProjectiveCoreTailScanCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveCoreTailScanCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
