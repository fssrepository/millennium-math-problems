#pragma once

#include "local_sld_projective_core_tail_ledger.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveCoreHeightScanRow {
    SpectralInteger maximum_height = 0;
    std::size_t fixed_family_cardinality = 0;
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

struct LocalSldProjectiveCoreHeightScanReport {
    int cutoff = 0;
    std::vector<LocalSldProjectiveCoreHeightScanRow> rows;
    SpectralReal fitted_open_signed_height_slope = 0.0L;
    SpectralReal fitted_open_absolute_height_slope = 0.0L;
    bool every_decomposition_exact = false;
    bool each_core_is_fixed_and_finite = false;
    bool uniform_height_tail_bound_proved = false;
    bool full_local_lemma_proved = false;
};

class LocalSldProjectiveCoreHeightScan {
public:
    [[nodiscard]] static LocalSldProjectiveCoreHeightScanReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        std::vector<SpectralInteger> maximum_heights,
        int threads = 12,
        bool exclude_signature_123 = false,
        bool exclude_triple_family = false);
};

struct LocalSldProjectiveCoreHeightScanCliOptions {
    std::string state_path;
    std::string certificate_path;
    std::vector<SpectralInteger> maximum_heights;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveCoreHeightScanCli {
public:
    [[nodiscard]] static LocalSldProjectiveCoreHeightScanCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveCoreHeightScanCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
