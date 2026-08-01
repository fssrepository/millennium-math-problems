#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>

namespace lemma {

struct HelicalAdversaryCliOptions {
    std::string state_path;
    std::string state_output_path;
    std::string certificate_path;
    std::string selection = "heterochiral";
    std::string mode = "trajectory";
    int cutoff = 0;
    int iterations = 8;
    int line_search_steps = 16;
    int trajectory_steps = 10;
    int restarts = 1;
    int workers = 12;
    SpectralReal initial_step = 0.1L;
    SpectralReal restart_mutation = 0.03L;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
    std::uint64_t seed = 20260801;
};

class HelicalAdversaryCli {
public:
    [[nodiscard]] static HelicalAdversaryCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
};

int run_helical_adversary(
    const HelicalAdversaryCliOptions& options, std::ostream& out);

}  // namespace lemma
