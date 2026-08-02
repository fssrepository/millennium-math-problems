#pragma once

#include "projective_triad_geometry.hpp"
#include "projective_core_family.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

using ProjectiveSignature = ProjectivePrimitiveSignature;

struct FiniteProjectiveFamilyClosureReport {
    int maximum_cutoff = 0;
    std::vector<ProjectiveTriadGeometryCertificate> members;
    std::size_t family_cardinality = 0;
    std::size_t total_ordered_role_count = 0;
    SpectralReal maximum_summed_input_degree_ratio = 0.0L;
    SpectralReal maximum_summed_target_degree_ratio = 0.0L;
    Rational union_incidence_degree_power{1};
    Rational bilinear_l2_frequency_power{3, 2};
    Rational internal_quartet_frequency_power{5};
    Rational target_frequency_power{11, 2};
    Rational frequency_gain{-1, 2};
    bool nonempty_fixed_finite_family = false;
    bool unique_primitive_signatures = false;
    bool every_member_triangle_feasible = false;
    bool union_input_degree_bound_proved = false;
    bool union_target_degree_bound_proved = false;
    bool union_bilinear_shell_bound_proved = false;
    bool complete_internal_self_cross_decomposition = false;
    bool internal_normalization_terms_bound_proved = false;
    bool cutoff_independent_internal_family_bound = false;
    bool core_tail_coupling_bound_proved = false;
    bool growing_tail_internal_bound_proved = false;
    bool uniform_over_growing_families_proved = false;
    bool full_local_lemma_proved = false;
};

class FiniteProjectiveFamilyClosure {
public:
    [[nodiscard]] static FiniteProjectiveFamilyClosureReport certify(
        int maximum_cutoff,
        std::vector<ProjectiveSignature> primitive_squared_lengths);
};

struct FiniteProjectiveFamilyClosureOptions {
    int maximum_cutoff = 8;
    SpectralInteger maximum_height = 0;
    std::vector<ProjectiveSignature> primitive_squared_lengths;
    std::string certificate_path;
};

class FiniteProjectiveFamilyClosureCli {
public:
    [[nodiscard]] static FiniteProjectiveFamilyClosureOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const FiniteProjectiveFamilyClosureOptions& options,
        std::ostream& out);
};

}  // namespace lemma
