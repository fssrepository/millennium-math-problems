#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightTransferRow {
    SpectralInteger core_maximum_height = 0;
    int cutoff = 0;
    SpectralReal open_power_one = 0.0L;
    SpectralReal core_tail_power_one = 0.0L;
    SpectralReal tail_internal_power_one = 0.0L;
    SpectralReal next_shell_diagonal_power_one = 0.0L;
    SpectralReal next_shell_fraction_of_open_absolute_sum = 0.0L;
    SpectralReal open_without_next_shell_diagonal = 0.0L;
    SpectralReal open_effective_height_pairs = 0.0L;
    SpectralReal open_signed_alignment = 0.0L;
    SpectralReal reconstruction_error = 0.0L;
    std::string state_path;
};

struct LocalSldProjectiveHeightTransferReport {
    std::vector<LocalSldProjectiveHeightTransferRow> rows;
    SpectralReal fitted_open_height_slope = 0.0L;
    SpectralReal fitted_next_shell_diagonal_height_slope = 0.0L;
    SpectralReal minimum_next_shell_fraction = 0.0L;
    SpectralReal maximum_next_shell_fraction = 0.0L;
    bool every_matrix_exact = false;
    bool finite_scale_transfer_is_not_a_proof = true;
    bool uniform_weighted_height_bound_proved = false;
};

struct LocalSldProjectiveHeightTransferCliOptions {
    std::vector<SpectralInteger> core_maximum_heights;
    std::vector<std::string> state_paths;
    std::string certificate_path;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveHeightTransferCli {
public:
    [[nodiscard]] static LocalSldProjectiveHeightTransferCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightTransferCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
