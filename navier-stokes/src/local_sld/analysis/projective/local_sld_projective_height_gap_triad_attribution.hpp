#pragma once

#include "spectral_state.hpp"
#include "triad_partition.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightGapTriadContribution {
    std::array<SpectralInteger, 3> primitive_shape{};
    WaveVector advecting_wave;
    WaveVector advected_wave;
    SpectralReal advecting_energy = 0.0L;
    SpectralReal advected_energy = 0.0L;
    SpectralReal coefficient_magnitude = 0.0L;
    SpectralReal contribution_norm2 = 0.0L;
    SpectralReal output_pairing = 0.0L;
    SpectralReal signed_output_fraction = 0.0L;
    SpectralReal absolute_output_fraction = 0.0L;
};

struct LocalSldProjectiveHeightGapTriadShellReport {
    int shell = 0;
    SpectralInteger minimum_height = 0;
    SpectralInteger maximum_height = 0;
    std::size_t shape_count = 0;
    std::size_t shell_interaction_count = 0;
    std::size_t target_interaction_count = 0;
    SpectralReal output_norm2 = 0.0L;
    SpectralReal reconstructed_output_norm2 = 0.0L;
    SpectralReal output_reconstruction_error = 0.0L;
    SpectralReal signed_output_fraction_sum = 0.0L;
    SpectralReal absolute_output_fraction_sum = 0.0L;
    SpectralReal effective_target_interaction_count = 0.0L;
    SpectralReal top_absolute_fraction_share = 0.0L;
    bool exact_output_reconstruction = false;
    std::vector<LocalSldProjectiveHeightGapTriadContribution>
        top_contributions;
};

struct LocalSldProjectiveHeightGapTriadAttributionReport {
    int cutoff = 0;
    WaveVector output_wave;
    LocalSldProjectiveHeightGapTriadShellReport first;
    LocalSldProjectiveHeightGapTriadShellReport second;
    bool finite = false;
};

// Attributes one shared Fourier output of two primitive-height shells to the
// exact ordered advection interactions which generate it.
class LocalSldProjectiveHeightGapTriadAttribution {
public:
    [[nodiscard]] static
    LocalSldProjectiveHeightGapTriadAttributionReport analyze(
        const SpectralState& state,
        TriadSelection selection,
        int first_shell,
        int second_shell,
        WaveVector output_wave,
        std::size_t top_interactions = 24,
        int threads = 12);
};

struct LocalSldProjectiveHeightGapTriadAttributionCliOptions {
    std::string state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    int first_shell = 0;
    int second_shell = 1;
    WaveVector output_wave;
    std::size_t top_interactions = 24;
    int threads = 12;
};

class LocalSldProjectiveHeightGapTriadAttributionCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveHeightGapTriadAttributionCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightGapTriadAttributionCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
