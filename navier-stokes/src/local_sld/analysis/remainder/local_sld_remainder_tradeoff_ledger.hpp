#pragma once

#include "local_sld_block_objective.hpp"
#include "local_sld_remainder_double_square.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldRemainderTradeoffRow {
    int cutoff = 0;
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal remainder_bracket = 0.0L;
    SpectralReal signed_lqc3_ratio = 0.0L;
    SpectralReal bracket_constant_ratio = 0.0L;
    SpectralReal normalized_stretching = 0.0L;
    SpectralReal signed_linear_tradeoff = 0.0L;
    SpectralReal power_one_commutator = 0.0L;
    SpectralReal power_one_advecting = 0.0L;
    SpectralReal power_one_enstrophy_normalization = 0.0L;
    SpectralReal power_one_palinstrophy_normalization = 0.0L;
    SpectralReal power_one_absolute_component_sum = 0.0L;
    SpectralReal power_one_cancellation_fraction = 0.0L;
    SpectralReal power_one_component_reconstruction_error = 0.0L;
    SpectralReal linear_tradeoff = 0.0L;
    SpectralReal quadratic_tradeoff = 0.0L;
    SpectralReal cubic_tradeoff = 0.0L;
    SpectralReal cubic_shape_numerator = 0.0L;
    SpectralReal shape_denominator = 0.0L;
    SpectralReal exact_shape_factor = 0.0L;
    SpectralReal linear_shape_multiplier = 0.0L;
    SpectralReal linear_reconstructed_block_ratio = 0.0L;
    SpectralReal factorized_block_ratio = 0.0L;
    SpectralReal direct_block_ratio = 0.0L;
    SpectralReal lqc3_to_block_depletion = 0.0L;
    SpectralReal upper_envelope_ratio = 0.0L;
    SpectralReal negative_square_ratio = 0.0L;
    SpectralReal shape_reconstruction_error = 0.0L;
    SpectralReal linear_reconstruction_error = 0.0L;
    SpectralReal product_reconstruction_error = 0.0L;
    bool scalar_multiplier_bound_satisfied = false;
    bool exact_factorization = false;
    bool cutoff_independent_tradeoff_proved = false;
};

struct LocalSldRemainderTradeoffReport {
    std::vector<LocalSldRemainderTradeoffRow> rows;
    SpectralReal bracket_constant_cutoff_slope = 0.0L;
    SpectralReal normalized_stretching_cutoff_slope = 0.0L;
    SpectralReal shape_factor_cutoff_slope = 0.0L;
    SpectralReal block_ratio_cutoff_slope = 0.0L;
    bool all_factorizations_exact = false;
    bool cutoff_independent_tradeoff_proved = false;
};

class LocalSldRemainderTradeoffLedger {
public:
    [[nodiscard]] static LocalSldRemainderTradeoffRow analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state);

    [[nodiscard]] static LocalSldRemainderTradeoffReport analyze(
        const SpectralDynamics& dynamics,
        const std::vector<SpectralState>& states);
};

struct LocalSldRemainderTradeoffCliOptions {
    std::vector<std::string> state_paths;
    std::string certificate_path;
    int threads = 12;
};

class LocalSldRemainderTradeoffCli {
public:
    [[nodiscard]] static LocalSldRemainderTradeoffCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldRemainderTradeoffCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
