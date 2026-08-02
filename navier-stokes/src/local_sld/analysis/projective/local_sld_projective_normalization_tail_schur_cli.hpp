#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationTailSchurCliOptions {
    std::vector<SpectralInteger> core_maximum_heights;
    std::vector<std::string> state_paths;
    std::string selection = "double-triple-remainder-without-123";
    std::string certificate_path;
    int threads = 12;
};

class LocalSldProjectiveNormalizationTailSchurCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveNormalizationTailSchurCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationTailSchurCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
