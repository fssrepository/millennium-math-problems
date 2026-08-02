#pragma once

#include "spectral_dynamics.hpp"
#include "triad_partition.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightGapOutputMode {
    WaveVector wave;
    SpectralReal first_h2_norm2 = 0.0L;
    SpectralReal second_h2_norm2 = 0.0L;
    SpectralReal pairing_real = 0.0L;
    SpectralReal pairing_imaginary = 0.0L;
    SpectralReal pairing_magnitude = 0.0L;
    SpectralReal modal_correlation = 0.0L;
    SpectralReal first_energy_fraction = 0.0L;
    SpectralReal second_energy_fraction = 0.0L;
    SpectralReal signed_gram_fraction = 0.0L;
    SpectralReal absolute_pairing_fraction = 0.0L;
};

struct LocalSldProjectiveHeightGapOutputReport {
    int cutoff = 0;
    int first_shell = 0;
    int second_shell = 0;
    int shell_gap = 0;
    std::size_t first_shape_count = 0;
    std::size_t second_shape_count = 0;
    std::size_t first_interaction_count = 0;
    std::size_t second_interaction_count = 0;
    std::size_t shared_conjugate_pair_count = 0;
    SpectralReal first_h2_norm2 = 0.0L;
    SpectralReal second_h2_norm2 = 0.0L;
    SpectralReal gram_pairing = 0.0L;
    SpectralReal reconstructed_gram_pairing = 0.0L;
    SpectralReal absolute_modal_pairing = 0.0L;
    SpectralReal correlation = 0.0L;
    SpectralReal first_shared_energy_fraction = 0.0L;
    SpectralReal second_shared_energy_fraction = 0.0L;
    SpectralReal effective_shared_mode_count = 0.0L;
    SpectralReal top_signed_gram_fraction = 0.0L;
    SpectralReal top_absolute_pairing_fraction = 0.0L;
    SpectralReal gram_reconstruction_error = 0.0L;
    bool exact_gram_reconstruction = false;
    bool finite = false;
    std::vector<LocalSldProjectiveHeightGapOutputMode> top_modes;
};

// Resolves the H2 Gram pairing between two primitive-height shells into
// conjugate Fourier-output pairs. This is a diagnostic ledger, not a bound.
class LocalSldProjectiveHeightGapOutputLedger {
public:
    [[nodiscard]] static LocalSldProjectiveHeightGapOutputReport analyze(
        const SpectralState& state,
        TriadSelection selection,
        int first_shell,
        int second_shell,
        std::size_t top_modes = 24,
        int threads = 12);
};

struct LocalSldProjectiveHeightGapOutputCliOptions {
    std::string state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    int first_shell = 0;
    int second_shell = 1;
    std::size_t top_modes = 24;
    int threads = 12;
};

class LocalSldProjectiveHeightGapOutputCli {
public:
    [[nodiscard]] static LocalSldProjectiveHeightGapOutputCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightGapOutputCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
