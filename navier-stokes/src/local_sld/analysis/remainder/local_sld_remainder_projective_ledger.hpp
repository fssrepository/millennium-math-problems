#pragma once

#include "local_sld_remainder_signature_ledger.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldRemainderProjectiveEntry {
    std::array<SpectralInteger, 3> primitive_squared_lengths{};
    std::size_t scale_count = 0;
    std::size_t interactions = 0;
    SpectralInteger minimum_scale = 0;
    SpectralInteger maximum_scale = 0;
    SpectralReal bracket_total = 0.0L;
    SpectralReal power_one_total = 0.0L;
    SpectralReal absolute_scale_sum = 0.0L;
    SpectralReal squared_scale_sum = 0.0L;
    SpectralReal effective_scales = 0.0L;
    SpectralReal signed_scale_alignment = 0.0L;
    SpectralReal absolute_fraction = 0.0L;
};

struct LocalSldRemainderProjectiveReport {
    int cutoff = 0;
    std::size_t signature_count = 0;
    std::size_t projective_shape_count = 0;
    bool excludes_signature_123 = false;
    bool excludes_triple_family = false;
    SpectralReal expected_power_one_total = 0.0L;
    SpectralReal reconstructed_power_one_total = 0.0L;
    SpectralReal reconstruction_error = 0.0L;
    SpectralReal absolute_projective_sum = 0.0L;
    SpectralReal squared_projective_sum = 0.0L;
    SpectralReal effective_projective_shapes = 0.0L;
    SpectralReal dominant_projective_fraction = 0.0L;
    SpectralReal signed_projective_alignment = 0.0L;
    bool exact_reconstruction = false;
    bool cutoff_independent_projective_sum_proved = false;
    std::vector<LocalSldRemainderProjectiveEntry> shapes;
};

class LocalSldRemainderProjectiveLedger {
public:
    [[nodiscard]] static LocalSldRemainderProjectiveReport analyze(
        const LocalSldRemainderSignatureReport& signatures);
};

struct LocalSldRemainderProjectiveCliOptions {
    std::string state_path;
    std::string certificate_path;
    int top = 64;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldRemainderProjectiveCli {
public:
    [[nodiscard]] static LocalSldRemainderProjectiveCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldRemainderProjectiveCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
