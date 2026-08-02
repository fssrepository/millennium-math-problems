#pragma once

#include "local_quartic_closure_objective.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldRemainderSignatureEntry {
    std::array<SpectralInteger, 3> squared_lengths{};
    std::size_t interactions = 0;
    SpectralReal stretching = 0.0L;
    SpectralReal palinstrophy_cross = 0.0L;
    SpectralReal outer_square = 0.0L;
    SpectralReal advected_commutator = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal palinstrophy_normalization = 0.0L;
    SpectralReal advecting_nested = 0.0L;
    SpectralReal total = 0.0L;
    SpectralReal absolute_fraction = 0.0L;
    SpectralReal target_ratio = 0.0L;
    SpectralReal power_one_ratio = 0.0L;
};

struct LocalSldRemainderSignatureReport {
    int cutoff = 0;
    std::size_t selected_interactions = 0;
    std::size_t signature_count = 0;
    int threads = 1;
    bool excludes_signature_123 = false;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal target_scale = 0.0L;
    SpectralReal full_bracket = 0.0L;
    SpectralReal full_target_ratio = 0.0L;
    SpectralReal common_full_stretching = 0.0L;
    SpectralReal power_one_scale = 0.0L;
    SpectralReal power_one_total = 0.0L;
    SpectralReal power_one_reconstruction_error = 0.0L;
    SpectralReal reconstructed_stretching = 0.0L;
    SpectralReal reconstructed_palinstrophy_cross = 0.0L;
    SpectralReal reconstructed_outer_square = 0.0L;
    SpectralReal reconstructed_advected_commutator = 0.0L;
    SpectralReal reconstructed_enstrophy_normalization = 0.0L;
    SpectralReal reconstructed_palinstrophy_normalization = 0.0L;
    SpectralReal reconstructed_advecting_nested = 0.0L;
    SpectralReal reconstructed_bracket = 0.0L;
    SpectralReal stretching_reconstruction_error = 0.0L;
    SpectralReal palinstrophy_cross_reconstruction_error = 0.0L;
    SpectralReal bracket_reconstruction_error = 0.0L;
    SpectralReal absolute_contribution_sum = 0.0L;
    SpectralReal squared_contribution_sum = 0.0L;
    SpectralReal effective_contributing_signatures = 0.0L;
    SpectralReal dominant_absolute_fraction = 0.0L;
    SpectralReal signed_cancellation_ratio = 0.0L;
    SpectralReal signed_amplification = 0.0L;
    bool exact_reconstruction = false;
    bool cutoff_independent_bound_proved = false;
    std::vector<LocalSldRemainderSignatureEntry> signatures;
};

class LocalSldRemainderSignatureLedger {
public:
    [[nodiscard]] static LocalSldRemainderSignatureReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        int threads = 12,
        bool exclude_signature_123 = false);
};

struct LocalSldRemainderSignatureCliOptions {
    std::string state_path;
    std::string certificate_path;
    int top = 64;
    int threads = 12;
    bool exclude_signature_123 = false;
};

class LocalSldRemainderSignatureCli {
public:
    [[nodiscard]] static LocalSldRemainderSignatureCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldRemainderSignatureCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
