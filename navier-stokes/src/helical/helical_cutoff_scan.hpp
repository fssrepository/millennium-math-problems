#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>

namespace lemma {

struct HelicalCutoffScanOptions {
    std::string state_path;
    std::string certificate_path;
    std::string selection = "heterochiral";
    std::string spread = "all";
    int minimum_cutoff = 3;
    int maximum_cutoff = 5;
    int trajectory_steps = 10;
    int workers = 12;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
    SpectralReal convergence_tolerance = 1.0e-6L;
    std::uint64_t seed = 20260801;
};

class HelicalCutoffScan {
public:
    [[nodiscard]] static HelicalCutoffScanOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const HelicalCutoffScanOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
