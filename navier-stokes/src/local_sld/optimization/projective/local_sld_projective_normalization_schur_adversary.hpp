#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldProjectiveNormalizationSchurAdversaryOptions {
    std::string state_path;
    std::string output_state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    SpectralInteger core_maximum_height = 8;
    int row_shell = 5;
    int iterations = 16;
    int line_search_steps = 14;
    int lbfgs_history = 8;
    int threads = 12;
    SpectralReal initial_step = 0.1L;
    std::uint64_t validation_seed = 20260802;
};

class LocalSldProjectiveNormalizationSchurAdversaryCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveNormalizationSchurAdversaryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationSchurAdversaryOptions& options,
        std::ostream& out);
};

}  // namespace lemma
