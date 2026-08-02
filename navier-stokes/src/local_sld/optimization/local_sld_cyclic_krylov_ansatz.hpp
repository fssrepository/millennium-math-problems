#pragma once

#include "local_sld_trajectory_adjoint.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldCyclicKrylovOptions {
    int cutoff = 3;
    int warm_angle_samples = 96;
    int iterations = 20;
    int line_search_steps = 14;
    int trajectory_steps = 500;
    int threads = 12;
    SpectralReal initial_step = 0.2L;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
    std::string backend = "auto";
    std::string warm_state_path;
    std::string certificate_path;
    std::string state_path;
};

struct LocalSldCyclicKrylovReport {
    SpectralState state;
    LocalSldTrajectoryValue value;
    LocalSldTrajectoryValue refined_value;
    std::array<SpectralReal, 3> coefficients{};
    std::array<SpectralReal, 3> energy_fractions{};
    SpectralReal maximum_gram_error = 0.0L;
    SpectralReal warm_projection_energy = 0.0L;
    SpectralReal warm_projection_residual = 0.0L;
    SpectralReal restricted_gradient_norm = 0.0L;
    SpectralReal projected_full_gradient_norm = 0.0L;
    SpectralReal time_step_relative_error = 0.0L;
    int accepted_steps = 0;
    int evaluations = 0;
    bool finite_search_is_not_a_proof = true;
};

class LocalSldCyclicKrylovAnsatz {
public:
    [[nodiscard]] static LocalSldCyclicKrylovReport optimize(
        const LocalSldCyclicKrylovOptions& options);
};

class LocalSldCyclicKrylovCli {
public:
    [[nodiscard]] static LocalSldCyclicKrylovOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldCyclicKrylovOptions& options,
        std::ostream& out);
};

}  // namespace lemma
