#pragma once

#include "spectral_dynamics.hpp"
#include "triad_partition.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationTailSchurShell {
    int shell = 0;
    SpectralInteger minimum_primitive_height = 0;
    SpectralInteger maximum_primitive_height = 0;
    std::size_t shape_count = 0;
    std::size_t interaction_count = 0;
    SpectralReal aggregate_h2_norm2 = 0.0L;
    SpectralReal diagonal_tail_fraction = 0.0L;
    SpectralReal normalized_absolute_gram_row_sum = 0.0L;
};

struct LocalSldProjectiveNormalizationTailSchurPair {
    int first_shell = 0;
    int second_shell = 0;
    int shell_gap = 0;
    SpectralReal gram_pairing = 0.0L;
    SpectralReal normalized_absolute_gram = 0.0L;
    SpectralReal half_decay_weighted_correlation = 0.0L;
};

struct LocalSldProjectiveNormalizationTailSchurGap {
    int shell_gap = 0;
    std::size_t pair_count = 0;
    SpectralReal maximum_normalized_absolute_gram = 0.0L;
    SpectralReal sum_normalized_absolute_gram = 0.0L;
    SpectralReal maximum_half_decay_weighted_correlation = 0.0L;
};

struct LocalSldProjectiveNormalizationTailSchurReport {
    int cutoff = 0;
    SpectralInteger core_maximum_height = 0;
    std::size_t selected_shape_count = 0;
    std::size_t tail_shape_count = 0;
    SpectralReal aggregate_tail_h2_norm2 = 0.0L;
    SpectralReal reconstructed_tail_h2_norm2 = 0.0L;
    SpectralReal diagonal_tail_h2_norm2 = 0.0L;
    SpectralReal absolute_gram_tail_h2 = 0.0L;
    SpectralReal maximum_normalized_absolute_gram_row_sum = 0.0L;
    SpectralReal schur_tail_h2_upper_bound = 0.0L;
    SpectralReal maximum_half_decay_weighted_correlation = 0.0L;
    SpectralReal half_decay_implied_row_bound = 0.0L;
    SpectralReal half_decay_tail_h2_upper_bound = 0.0L;
    SpectralReal aggregate_to_diagonal_ratio = 0.0L;
    SpectralReal aggregate_to_absolute_gram_ratio = 0.0L;
    SpectralReal aggregate_to_schur_bound_ratio = 0.0L;
    SpectralReal tail_h2_reconstruction_error = 0.0L;
    SpectralReal objective_tail_h2_reconstruction_error = 0.0L;
    SpectralReal full_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal selected_aggregate_h1_norm2 = 0.0L;
    SpectralReal normalization_common_factor = 0.0L;
    SpectralReal squared_cauchy_majorant = 0.0L;
    SpectralReal diagonal_squared_majorant = 0.0L;
    SpectralReal schur_squared_majorant = 0.0L;
    SpectralReal height_half_compensated_squared_cauchy_majorant = 0.0L;
    SpectralReal height_half_compensated_diagonal_squared_majorant = 0.0L;
    SpectralReal height_half_compensated_schur_squared_majorant = 0.0L;
    SpectralReal height_half_compensated_half_decay_squared_majorant = 0.0L;
    bool exact_tail_h2_reconstruction = false;
    bool finite_schur_inequality_verified = false;
    bool finite_half_gap_decay_inequality_verified = false;
    bool finite = false;
    std::vector<LocalSldProjectiveNormalizationTailSchurShell> shells;
    std::vector<LocalSldProjectiveNormalizationTailSchurPair> pairs;
    std::vector<LocalSldProjectiveNormalizationTailSchurGap> gaps;
};

// Exact finite Gram/Schur factorization of
// T2 = ||A sum_(primitive height > H) B_shape(u,u)||_2^2.
// It separates the diagonal height-shell mass from inter-shell coherence.
class LocalSldProjectiveNormalizationTailSchurLedger {
public:
    [[nodiscard]] static
    LocalSldProjectiveNormalizationTailSchurReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        int threads = 12);
};

}  // namespace lemma
