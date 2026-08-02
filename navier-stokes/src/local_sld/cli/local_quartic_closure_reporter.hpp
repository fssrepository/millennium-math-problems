#pragma once

#include "local_quartic_closure_adversary.hpp"

#include <iosfwd>

namespace lemma {

class LocalQuarticClosureReporter {
public:
    static void write_artifacts(
        LocalQuarticClosureAdversaryReport& report,
        const LocalQuarticClosureAdversaryOptions& options);
    static void print_summary(
        const LocalQuarticClosureAdversaryReport& report,
        std::ostream& out);
};

}  // namespace lemma
