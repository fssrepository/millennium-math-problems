#pragma once

#include "proof_scaling.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct ProjectiveSquareFunctionClosureReport {
    Rational uniform_ray_incidence_degree_power{1};
    Rational derivative_frequency_power{1};
    Rational squared_function_squared_frequency_power{3};
    Rational squared_function_frequency_power{3, 2};
    Rational candidate_diagonal_quartet_frequency_power{5};
    Rational target_frequency_power{11, 2};
    Rational candidate_diagonal_frequency_gain{-1, 2};
    bool primitive_ray_partition_is_disjoint = false;
    bool uniform_single_ray_incidence_constant_proved = false;
    bool bilinear_projective_square_function_bound_proved = false;
    bool square_function_has_target_power_gain = false;
    bool diagonal_projective_quartet_sum_proved = false;
    bool coherent_projective_synthesis_bound_proved = false;
    bool cutoff_independent_projective_synthesis_bound_rejected = false;
    bool coherent_fan_zero_stretching_proved = false;
    bool cross_ray_quartet_bound_proved = false;
    bool power_one_tradeoff_bound_proved = false;
    bool full_local_lemma_proved = false;
};

class ProjectiveSquareFunctionClosure {
public:
    [[nodiscard]] static ProjectiveSquareFunctionClosureReport certify();
};

struct ProjectiveSquareFunctionClosureOptions {
    std::string certificate_path;
};

class ProjectiveSquareFunctionClosureCli {
public:
    [[nodiscard]] static ProjectiveSquareFunctionClosureOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const ProjectiveSquareFunctionClosureOptions& options,
        std::ostream& out);
};

}  // namespace lemma
