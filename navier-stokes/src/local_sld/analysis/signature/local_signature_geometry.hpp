#pragma once

#include "proof_scaling.hpp"
#include "spectral_state.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSignatureCountRow {
    int cutoff = 0;
    std::size_t ordered_signatures = 0;
    std::size_t maximum_input_degree = 0;
    std::size_t maximum_target_degree = 0;
    std::size_t ordered_signature_degree_bound = 0;
};

struct LocalSignatureGeometryCertificate {
    int maximum_cutoff = 0;
    SpectralReal maximum_input_degree_ratio = 0.0L;
    SpectralReal maximum_target_degree_ratio = 0.0L;
    bool all_fixed_signature_degree_bounds_hold = true;
    std::vector<LocalSignatureCountRow> rows;
};

struct SignatureFamilyClosureRow {
    Rational signature_count_power{0};
    Rational union_degree_power{1};
    Rational transfer_frequency_power{7, 2};
    Rational transfer_to_viscosity_frequency_power{-1, 2};
    bool energy_level_high_frequency_absorption = true;
};

struct SignatureFamilyClosure {
    Rational derivative_weight_power{3};
    Rational fixed_signature_degree_power{1};
    Rational cauchy_degree_power{1, 2};
    Rational viscous_frequency_power{4};
    Rational critical_signature_count_power{1};
    Rational critical_signed_amplification_power{1, 2};
    SignatureFamilyClosureRow finite_signature_family;
    SignatureFamilyClosureRow critical_signature_family;
    SignatureFamilyClosureRow dense_signature_family;
    bool square_summed_fixed_signature_bound = true;
    bool effective_count_replaces_raw_count = true;
    bool signed_amplification_preserves_cancellation = true;
    bool closing_requires_sublinear_signature_count = true;
};

class LocalSignatureGeometry {
public:
    [[nodiscard]] static LocalSignatureGeometryCertificate certify(
        int maximum_cutoff);
    [[nodiscard]] static SignatureFamilyClosure analyze_closure();
};

struct LocalSignatureCliOptions {
    int maximum_cutoff = 5;
    std::string certificate_path;
};

class LocalSignatureCli {
public:
    [[nodiscard]] static LocalSignatureCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const LocalSignatureCliOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
