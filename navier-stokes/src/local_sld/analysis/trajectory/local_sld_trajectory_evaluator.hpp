#pragma once

#include "local_sld_trajectory_adjoint.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldTrajectoryEvaluatorOptions {
    std::string state_path;
    std::string certificate_path;
    std::string backend = "auto";
    int trajectory_steps = 500;
    int threads = 12;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
};

struct LocalSldTrajectoryEvaluatorReport {
    LocalSldTrajectoryValue initial;
    LocalSldTrajectoryValue terminal;
    LocalSldTrajectoryValue maximum;
    LocalSldTrajectoryValue refined_maximum;
    SpectralReal time_step_relative_error = 0.0L;
    int cutoff = 0;
    bool finite = false;
    bool finite_evaluation_is_not_a_proof = true;
};

class LocalSldTrajectoryEvaluator {
public:
    [[nodiscard]] static LocalSldTrajectoryEvaluatorReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        SpectralReal viscosity,
        SpectralReal time_step,
        int trajectory_steps);
};

class LocalSldTrajectoryEvaluatorCli {
public:
    [[nodiscard]] static LocalSldTrajectoryEvaluatorOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldTrajectoryEvaluatorOptions& options,
        std::ostream& out);
};

}  // namespace lemma
