#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>

namespace lemma {

struct SpectralStateTransformOptions {
    std::string state_path;
    std::string output_path;
    int target_cutoff = 0;
    std::uint64_t seed = 20260802;
};

class SpectralStateTransform {
public:
    [[nodiscard]] static SpectralState to_cutoff(
        const SpectralState& state,
        int target_cutoff,
        std::uint64_t seed);
};

class SpectralStateTransformCli {
public:
    [[nodiscard]] static SpectralStateTransformOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const SpectralStateTransformOptions& options,
        std::ostream& out);
};

}  // namespace lemma
