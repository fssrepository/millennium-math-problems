#pragma once

#include "proof_scaling.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct RemainderQuartetClosureReport {
    Rational dense_incidence_degree_power{3};
    Rational dense_bilinear_frequency_power{5, 2};
    Rational dense_stretching_frequency_power{9, 2};
    Rational dense_weighted_stretching_frequency_power{13, 2};
    Rational dense_structural_bracket_frequency_power{7};
    Rational dense_normalization_bracket_frequency_power{7};
    Rational target_frequency_power{11, 2};
    Rational dense_frequency_loss{3, 2};
    Rational required_bilinear_frequency_power{7, 4};
    Rational required_effective_incidence_degree_power{3, 2};
    Rational required_incidence_reduction_power{3, 2};
    Rational fixed_signature_incidence_degree_power{1};
    Rational fixed_signature_bilinear_frequency_power{3, 2};
    Rational fixed_signature_bracket_frequency_power{5};
    Rational fixed_signature_frequency_gain{-1, 2};
    bool dense_energy_only_count_closes = false;
    bool every_fixed_signature_closes = false;
    bool remainder_requires_collective_cancellation = false;
    bool one_sided_double_square_reduction = false;
    bool stretching_vjp_commutator_identity = false;
    bool standalone_commutator_envelope_bound_proved = false;
    bool commutator_absorption_bound_proved = false;
    bool cutoff_independent_remainder_bound_proved = false;
    bool full_local_lemma_proved = false;
};

class RemainderQuartetClosure {
public:
    [[nodiscard]] static RemainderQuartetClosureReport certify();
};

struct RemainderQuartetClosureOptions {
    std::string certificate_path;
};

class RemainderQuartetClosureCli {
public:
    [[nodiscard]] static RemainderQuartetClosureOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const RemainderQuartetClosureOptions& options,
        std::ostream& out);
};

}  // namespace lemma
