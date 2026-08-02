#pragma once

#include "local_quartic_closure_adversary.hpp"

#include <iosfwd>

namespace lemma {

class LocalQuarticClosureCli {
public:
    [[nodiscard]] static LocalQuarticClosureAdversaryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalQuarticClosureAdversaryOptions& options,
        std::ostream& out);
};

}  // namespace lemma
