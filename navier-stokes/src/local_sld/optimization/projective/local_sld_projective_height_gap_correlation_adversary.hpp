#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldProjectiveHeightGapCorrelationAdversaryOptions {
    std::string state_path;
    std::string output_state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    int first_shell = 4;
    int second_shell = 8;
    int iterations = 24;
    int line_search_steps = 14;
    int lbfgs_history = 8;
    int threads = 12;
    SpectralReal initial_step = 0.1L;
    std::uint64_t validation_seed = 20260802;
};

class LocalSldProjectiveHeightGapCorrelationAdversaryCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveHeightGapCorrelationAdversaryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveHeightGapCorrelationAdversaryOptions& options,
        std::ostream& out);
};

}  // namespace lemma
