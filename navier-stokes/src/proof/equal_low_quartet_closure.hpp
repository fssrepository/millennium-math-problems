#pragma once

#include "equal_low_triad_geometry.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct EqualLowQuartetClosureReport {
    int target_length_multiplier = 3;
    EqualLowTriadGeometryCertificate geometry;
    Rational bilinear_l2_frequency_power{3, 2};
    Rational bracket_frequency_power{5};
    Rational target_frequency_power{11, 2};
    Rational frequency_gain{-1, 2};
    bool fixed_angle_plane_sphere_bound_proved = false;
    bool closed_single_shell_power_bound = false;
    bool structural_entries_neighbor_shell_local = false;
    bool direct_normalization_target_bound_proved = false;
    bool cutoff_independent_closed_family_bound = false;
    bool full_local_lemma_proved = false;
};

class EqualLowQuartetClosure {
public:
    [[nodiscard]] static EqualLowQuartetClosureReport certify(
        int maximum_cutoff,
        int target_length_multiplier = 3);
};

struct EqualLowQuartetClosureOptions {
    int maximum_cutoff = 8;
    int target_length_multiplier = 3;
    std::string certificate_path;
};

class EqualLowQuartetClosureCli {
public:
    [[nodiscard]] static EqualLowQuartetClosureOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const EqualLowQuartetClosureOptions& options,
        std::ostream& out);
};

}  // namespace lemma
