#pragma once

#include "local_quartic_closure_objective.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldCyclicAnsatzOptions {
    int coarse_samples = 4096;
    int refinement_iterations = 96;
    std::string certificate_path;
    std::string state_path;
};

struct LocalSldCyclicAnsatzReport {
    SpectralState state;
    LocalQuarticClosureObjectiveValue value;
    LocalQuarticClosureObjectiveValue pure_axis_value;
    LocalQuarticClosureObjectiveValue pure_response_value;
    SpectralReal angle = 0.0L;
    SpectralReal axis_energy_fraction = 0.0L;
    SpectralReal response_energy_fraction = 0.0L;
    SpectralReal basis_inner_product = 0.0L;
    SpectralReal projected_gradient_norm = 0.0L;
    SpectralReal pure_axis_identity_error = 0.0L;
    int coarse_samples = 0;
    int refinement_iterations = 0;
    bool finite_search_is_not_a_proof = true;
};

class LocalSldCyclicAnsatz {
public:
    [[nodiscard]] static LocalSldCyclicAnsatzReport optimize(
        const LocalSldCyclicAnsatzOptions& options);
};

class LocalSldCyclicAnsatzCli {
public:
    [[nodiscard]] static LocalSldCyclicAnsatzOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldCyclicAnsatzOptions& options,
        std::ostream& out);
};

}  // namespace lemma
