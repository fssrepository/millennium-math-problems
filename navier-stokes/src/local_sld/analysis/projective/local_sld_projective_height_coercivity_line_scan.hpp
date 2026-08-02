#pragma once

#include "spectral_dynamics.hpp"
#include "spectral_state.hpp"
#include "triad_partition.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightCoercivityLineRow {
    SpectralReal parameter = 0.0L;
    SpectralReal paired_envelope = 0.0L;
    SpectralReal outer_h1_weight = 0.0L;
    SpectralReal coercivity_ratio = 0.0L;
    SpectralReal squared_coercivity_ratio = 0.0L;
};

struct LocalSldProjectiveHeightCoercivityLineReport {
    int cutoff = 0;
    SpectralReal maximum_ratio = 0.0L;
    SpectralReal maximum_ratio_parameter = 0.0L;
    SpectralReal minimum_outer_h1_weight = 0.0L;
    SpectralReal minimum_outer_parameter = 0.0L;
    SpectralReal low_outer_paired_exponent = 0.0L;
    std::vector<LocalSldProjectiveHeightCoercivityLineRow> rows;
};

class LocalSldProjectiveHeightCoercivityLineScan {
public:
    [[nodiscard]] static LocalSldProjectiveHeightCoercivityLineReport
    analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& left,
        const SpectralState& right,
        TriadSelection selection,
        SpectralReal minimum_parameter,
        SpectralReal maximum_parameter,
        int samples,
        int threads);
};

struct LocalSldProjectiveHeightCoercivityLineScanOptions {
    std::string left_state_path;
    std::string right_state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder";
    SpectralReal minimum_parameter = 0.0L;
    SpectralReal maximum_parameter = 4.0L;
    int samples = 25;
    int threads = 12;
};

class LocalSldProjectiveHeightCoercivityLineScanCli {
public:
    [[nodiscard]] static LocalSldProjectiveHeightCoercivityLineScanOptions
    parse(int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightCoercivityLineScanOptions& options,
        std::ostream& out);
};

}  // namespace lemma
