#pragma once

#include "local_sld_projective_quartic_cross_ledger.hpp"
#include "local_sld_remainder_projective_ledger.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveCrossAttributionEntry {
    std::array<SpectralInteger, 3> primitive_squared_lengths{};
    SpectralReal full_attributed_power_one = 0.0L;
    SpectralReal self_power_one = 0.0L;
    SpectralReal cross_attributed_power_one = 0.0L;
    SpectralReal absolute_fraction = 0.0L;
};

struct LocalSldProjectiveCrossAttributionReport {
    int cutoff = 0;
    std::size_t projective_shape_count = 0;
    SpectralReal expected_cross_power_one = 0.0L;
    SpectralReal reconstructed_cross_power_one = 0.0L;
    SpectralReal reconstruction_error = 0.0L;
    SpectralReal absolute_cross_attribution_sum = 0.0L;
    SpectralReal squared_cross_attribution_sum = 0.0L;
    SpectralReal effective_cross_attribution_shapes = 0.0L;
    SpectralReal dominant_cross_attribution_fraction = 0.0L;
    SpectralReal signed_cross_attribution_alignment = 0.0L;
    bool exact_reconstruction = false;
    bool finite_ledger_is_not_a_proof = true;
    std::vector<LocalSldProjectiveCrossAttributionEntry> shapes;
};

class LocalSldProjectiveCrossAttribution {
public:
    [[nodiscard]] static LocalSldProjectiveCrossAttributionReport analyze(
        const LocalSldRemainderProjectiveReport& full_attribution,
        const LocalSldProjectiveQuarticCrossReport& self_cross);
};

struct LocalSldProjectiveCrossAttributionCliOptions {
    std::string state_path;
    std::string certificate_path;
    int top = 64;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveCrossAttributionCli {
public:
    [[nodiscard]] static LocalSldProjectiveCrossAttributionCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveCrossAttributionCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
