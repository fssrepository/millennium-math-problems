#pragma once

#include "projective_fan_geometry.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct ProjectiveFanObstructionOptions {
    int maximum_cutoff = 1024;
    std::string certificate_path;
};

class ProjectiveFanObstructionCli {
public:
    [[nodiscard]] static ProjectiveFanObstructionOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const ProjectiveFanObstructionOptions& options,
        std::ostream& out);
};

}  // namespace lemma
