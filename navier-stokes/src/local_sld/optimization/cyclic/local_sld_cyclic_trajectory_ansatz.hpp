#pragma once

#include "local_sld_trajectory_adjoint.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldCyclicTrajectoryOptions {
    int cutoff = 2;
    int coarse_samples = 256;
    int refinement_iterations = 64;
    int trajectory_steps = 500;
    int threads = 12;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
    std::string backend = "auto";
    std::string certificate_path;
    std::string state_path;
};

struct LocalSldCyclicTrajectoryReport {
    SpectralState state;
    LocalSldTrajectoryValue value;
    LocalSldTrajectoryValue refined_value;
    LocalQuarticClosureObjectiveValue initial_value;
    SpectralReal angle = 0.0L;
    SpectralReal axis_energy_fraction = 0.0L;
    SpectralReal response_energy_fraction = 0.0L;
    SpectralReal basis_inner_product = 0.0L;
    SpectralReal restricted_gradient = 0.0L;
    SpectralReal projected_full_gradient_norm = 0.0L;
    SpectralReal time_step_relative_error = 0.0L;
    int coarse_samples = 0;
    int refinement_iterations = 0;
    bool finite_search_is_not_a_proof = true;
};

class LocalSldCyclicTrajectoryAnsatz {
public:
    [[nodiscard]] static LocalSldCyclicTrajectoryReport optimize(
        const LocalSldCyclicTrajectoryOptions& options);
};

class LocalSldCyclicTrajectoryCli {
public:
    [[nodiscard]] static LocalSldCyclicTrajectoryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldCyclicTrajectoryOptions& options,
        std::ostream& out);
};

}  // namespace lemma
