#pragma once

#include "local_quartic_closure_objective.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveQuarticSelfEntry {
    std::array<SpectralInteger, 3> primitive_squared_lengths{};
    std::size_t interactions = 0;
    SpectralReal outer_square = 0.0L;
    SpectralReal advected_commutator = 0.0L;
    SpectralReal advecting_nested = 0.0L;
    SpectralReal enstrophy_normalization = 0.0L;
    SpectralReal palinstrophy_normalization = 0.0L;
    SpectralReal self_bracket = 0.0L;
    SpectralReal power_one_self = 0.0L;
    SpectralReal absolute_fraction = 0.0L;
};

struct LocalSldProjectiveQuarticCrossReport {
    int cutoff = 0;
    int threads = 1;
    bool excludes_signature_123 = false;
    bool excludes_triple_family = false;
    std::size_t projective_shape_count = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal full_bracket = 0.0L;
    SpectralReal full_outer_commutator = 0.0L;
    SpectralReal full_advecting_nested = 0.0L;
    SpectralReal full_enstrophy_normalization = 0.0L;
    SpectralReal full_palinstrophy_normalization = 0.0L;
    SpectralReal diagonal_outer_commutator = 0.0L;
    SpectralReal diagonal_advecting_nested = 0.0L;
    SpectralReal diagonal_enstrophy_normalization = 0.0L;
    SpectralReal diagonal_palinstrophy_normalization = 0.0L;
    SpectralReal diagonal_bracket = 0.0L;
    SpectralReal cross_bracket = 0.0L;
    SpectralReal power_one_scale = 0.0L;
    SpectralReal full_power_one = 0.0L;
    SpectralReal diagonal_power_one = 0.0L;
    SpectralReal cross_power_one = 0.0L;
    SpectralReal absolute_diagonal_power_one_sum = 0.0L;
    SpectralReal effective_diagonal_shapes = 0.0L;
    SpectralReal dominant_diagonal_fraction = 0.0L;
    SpectralReal diagonal_signed_alignment = 0.0L;
    SpectralReal full_component_reconstruction_error = 0.0L;
    SpectralReal bracket_decomposition_error = 0.0L;
    bool exact_decomposition = false;
    bool finite_ledger_is_not_a_proof = true;
    std::vector<LocalSldProjectiveQuarticSelfEntry> shapes;
};

class LocalSldProjectiveQuarticCrossLedger {
public:
    [[nodiscard]] static LocalSldProjectiveQuarticCrossReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        int threads = 12,
        bool exclude_signature_123 = false,
        bool exclude_triple_family = false);
};

struct LocalSldProjectiveQuarticCrossCliOptions {
    std::string state_path;
    std::string certificate_path;
    int top = 64;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveQuarticCrossCli {
public:
    [[nodiscard]] static LocalSldProjectiveQuarticCrossCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveQuarticCrossCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
