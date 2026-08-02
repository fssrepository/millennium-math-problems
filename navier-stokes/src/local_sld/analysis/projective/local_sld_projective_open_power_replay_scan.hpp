#pragma once

#include "local_sld_projective_core_tail_ledger.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveOpenPowerReplayRow {
    SpectralInteger core_maximum_height = 0;
    int cutoff = 0;
    std::size_t fixed_core_shape_count = 0;
    SpectralReal absolute_open_power_one = 0.0L;
    SpectralReal squared_open_power_one = 0.0L;
    SpectralReal core_tail_power_one = 0.0L;
    SpectralReal tail_internal_power_one = 0.0L;
    SpectralReal reconstruction_error = 0.0L;
    std::string state_path;
};

struct LocalSldProjectiveOpenPowerReplayReport {
    std::vector<LocalSldProjectiveOpenPowerReplayRow> rows;
    SpectralReal fitted_absolute_height_slope = 0.0L;
    SpectralReal fitted_squared_height_slope = 0.0L;
    bool every_reconstruction_exact = false;
    bool finite_optimized_replay_is_not_a_proof = true;
    bool uniform_height_tail_bound_proved = false;
    bool full_local_lemma_proved = false;
};

struct LocalSldProjectiveOpenPowerReplayCliOptions {
    std::vector<SpectralInteger> core_maximum_heights;
    std::vector<std::string> state_paths;
    std::string certificate_path;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveOpenPowerReplayCli {
public:
    [[nodiscard]] static LocalSldProjectiveOpenPowerReplayCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveOpenPowerReplayCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
