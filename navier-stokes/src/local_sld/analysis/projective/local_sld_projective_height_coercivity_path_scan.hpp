#pragma once

#include "spectral_dynamics.hpp"
#include "spectral_state.hpp"
#include "triad_partition.hpp"

#include <array>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightCoercivityPathRow {
    SpectralReal epsilon = 0.0L;
    SpectralReal energy = 0.0L;
    SpectralReal paired_envelope = 0.0L;
    SpectralReal outer_h1_weight = 0.0L;
    SpectralReal coercivity_ratio = 0.0L;
    SpectralReal squared_coercivity_ratio = 0.0L;
    SpectralReal signed_bracket = 0.0L;
    std::array<SpectralReal, 5> component_envelopes{};
};

struct LocalSldProjectiveHeightCoercivityPathReport {
    int cutoff = 0;
    WaveVector base_wave{};
    std::size_t shell_count = 0;
    SpectralReal base_selected_advection_h1_weight = 0.0L;
    SpectralReal paired_envelope_tail_slope = 0.0L;
    SpectralReal outer_h1_weight_tail_slope = 0.0L;
    SpectralReal ratio_tail_slope = 0.0L;
    bool base_is_selected_advection_null = false;
    std::vector<LocalSldProjectiveHeightCoercivityPathRow> rows;
};

class LocalSldProjectiveHeightCoercivityPathScan {
public:
    [[nodiscard]] static SpectralState monochromatic_base(
        const SpectralState& layout,
        WaveVector wave);

    [[nodiscard]] static LocalSldProjectiveHeightCoercivityPathReport
    analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& perturbation,
        TriadSelection selection,
        WaveVector base_wave,
        SpectralReal minimum_epsilon,
        SpectralReal maximum_epsilon,
        int samples,
        int threads);
};

struct LocalSldProjectiveHeightCoercivityPathScanOptions {
    std::string state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder";
    WaveVector base_wave{1, 0, 0};
    SpectralReal minimum_epsilon = 1e-6L;
    SpectralReal maximum_epsilon = 1e-1L;
    int samples = 13;
    int threads = 12;
};

class LocalSldProjectiveHeightCoercivityPathScanCli {
public:
    [[nodiscard]] static LocalSldProjectiveHeightCoercivityPathScanOptions
    parse(int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightCoercivityPathScanOptions& options,
        std::ostream& out);
};

}  // namespace lemma
