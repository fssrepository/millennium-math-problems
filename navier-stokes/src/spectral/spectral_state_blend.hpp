#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct SpectralStateBlendOptions {
    std::string left_state_path;
    std::string right_state_path;
    std::string output_path;
    SpectralReal right_weight = 0.01L;
};

class SpectralStateBlend {
public:
    [[nodiscard]] static SpectralState blend_on_energy_sphere(
        const SpectralState& left,
        const SpectralState& right,
        SpectralReal right_weight);

    [[nodiscard]] static SpectralState affine_normalized(
        const SpectralState& left,
        const SpectralState& right,
        SpectralReal parameter);
};

class SpectralStateBlendCli {
public:
    [[nodiscard]] static SpectralStateBlendOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const SpectralStateBlendOptions& options,
        std::ostream& out);
};

}  // namespace lemma
