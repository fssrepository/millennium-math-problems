#pragma once

#include "projective_triad_geometry.hpp"

#include <array>
#include <iosfwd>
#include <string>

namespace lemma {

struct ProjectiveQuartetClosureReport {
    ProjectiveTriadGeometryCertificate geometry;
    Rational incidence_degree_power{1};
    Rational bilinear_l2_frequency_power{3, 2};
    Rational bracket_frequency_power{5};
    Rational target_frequency_power{11, 2};
    Rational frequency_gain{-1, 2};
    bool fixed_projective_ray_incidence_proved = false;
    bool fixed_projective_ray_shell_bound_proved = false;
    bool fixed_projective_ray_normalization_bound_proved = false;
    bool cutoff_independent_fixed_projective_ray_bound = false;
    bool uniform_sum_over_projective_shapes_proved = false;
    bool full_local_lemma_proved = false;
};

class ProjectiveQuartetClosure {
public:
    [[nodiscard]] static ProjectiveQuartetClosureReport certify(
        int maximum_cutoff,
        std::array<SpectralInteger, 3> primitive_squared_lengths);
};

struct ProjectiveQuartetClosureOptions {
    int maximum_cutoff = 8;
    std::array<SpectralInteger, 3> primitive_squared_lengths{2, 3, 5};
    std::string certificate_path;
};

class ProjectiveQuartetClosureCli {
public:
    [[nodiscard]] static ProjectiveQuartetClosureOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const ProjectiveQuartetClosureOptions& options,
        std::ostream& out);
};

}  // namespace lemma
