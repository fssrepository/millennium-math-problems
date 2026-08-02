#pragma once

#include "orthogonal_triad_geometry.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct DoublingQuartetClosureReport {
    OrthogonalTriadGeometryCertificate geometry;
    Rational derivative_power{1};
    Rational incidence_degree_power{1};
    Rational bilinear_l2_frequency_power{3, 2};
    Rational bilinear_energy_power{1};
    Rational stretching_frequency_power{7, 2};
    Rational bracket_frequency_power{5};
    Rational bracket_energy_power{2};
    Rational required_frequency_power{11, 2};
    Rational required_energy_power{2};
    Rational frequency_gain{-1, 2};
    Rational global_enstrophy_power{3, 2};
    Rational global_palinstrophy_power{1, 2};
    Rational stretching_enstrophy_power{5, 4};
    Rational stretching_palinstrophy_power{1, 4};
    Rational weighted_stretching_enstrophy_power{1, 4};
    Rational weighted_stretching_palinstrophy_power{5, 4};
    Rational normalization_enstrophy_power{3, 2};
    Rational normalization_palinstrophy_power{1, 2};
    Rational target_enstrophy_power{5, 4};
    Rational target_palinstrophy_power{3, 4};
    bool target_energy_homogeneity_matches = false;
    bool every_quartet_entry_has_same_power = false;
    bool shell_sum_closes_by_zp = false;
    bool torus_spectral_gap_closes_target = false;
    bool closed_single_shell_power_bound = false;
    bool structural_entries_neighbor_shell_local = false;
    bool structural_entry_global_bound_proved = false;
    bool projected_normalization_bound_proved = false;
    bool direct_normalization_target_bound_proved = false;
    SpectralReal maximum_tested_normalization_ratio = 0.0L;
    int maximizing_first_shell = 0;
    int maximizing_second_shell = 0;
    SpectralReal maximizing_energy_ratio = 0.0L;
    bool normalization_sequence_screen_survives = false;
    bool normalization_sequence_bound_proved = false;
    Rational two_scale_energy_decay_power{11, 4};
    Rational two_scale_stretching_high_power{-5, 8};
    Rational two_scale_weighted_stretching_high_power{11, 8};
    Rational two_scale_palinstrophy_high_power{5, 4};
    Rational two_scale_shell_quartic_high_power{-1, 2};
    Rational two_scale_ratio_growth_power{1, 8};
    bool naive_cross_shell_bound_rejected = false;
    bool cutoff_independent_closed_family_bound = false;
    bool full_local_lemma_proved = false;
};

class DoublingQuartetClosure {
public:
    [[nodiscard]] static DoublingQuartetClosureReport certify(
        int maximum_cutoff);
};

struct DoublingQuartetClosureOptions {
    int maximum_cutoff = 8;
    std::string certificate_path;
};

class DoublingQuartetClosureCli {
public:
    [[nodiscard]] static DoublingQuartetClosureOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const DoublingQuartetClosureOptions& options,
        std::ostream& out);
};

}  // namespace lemma
