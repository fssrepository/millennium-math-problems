#pragma once

#include "local_sld_projective_height_matrix.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightSchurEntry {
    int first_shell = 0;
    int second_shell = 0;
    int shell_gap = 0;
    SpectralReal power_one = 0.0L;
    SpectralReal absolute_component_power_one_envelope = 0.0L;
    SpectralReal geometric_diagonal_scale = 0.0L;
    SpectralReal symmetric_geometric_ratio = 0.0L;
    bool has_nonzero_diagonal_scale = false;
};

struct LocalSldProjectiveHeightSchurGapRow {
    int shell_gap = 0;
    std::size_t pair_count = 0;
    std::size_t unscaled_pair_count = 0;
    SpectralReal signed_power_one_sum = 0.0L;
    SpectralReal absolute_power_one_sum = 0.0L;
    SpectralReal absolute_component_envelope_sum = 0.0L;
    SpectralReal maximum_symmetric_geometric_ratio = 0.0L;
    SpectralReal commutator_paired_envelope_sum = 0.0L;
    SpectralReal
        commutator_paired_outer_maximum_symmetric_geometric_ratio = 0.0L;
    SpectralReal commutator_term_envelope_sum = 0.0L;
    SpectralReal remainder_terms_envelope_sum = 0.0L;
    SpectralReal commutator_outer_maximum_symmetric_geometric_ratio = 0.0L;
    SpectralReal remainder_outer_maximum_symmetric_geometric_ratio = 0.0L;
    std::size_t commutator_paired_outer_unscaled_pair_count = 0;
};

struct LocalSldProjectiveHeightComponentSchurRow {
    const char* component = "";
    SpectralReal total_power_one_envelope = 0.0L;
    SpectralReal diagonal_power_one_envelope = 0.0L;
    SpectralReal maximum_off_diagonal_geometric_ratio = 0.0L;
    SpectralReal maximum_weighted_row_sum = 0.0L;
    std::size_t unscaled_off_diagonal_pair_count = 0;
    SpectralReal maximum_outer_weighted_diagonal_ratio = 0.0L;
    SpectralReal maximum_outer_weighted_off_diagonal_ratio = 0.0L;
    SpectralReal maximum_outer_weighted_row_sum = 0.0L;
    std::size_t outer_unscaled_pair_count = 0;
};

struct LocalSldProjectiveHeightSchurReport {
    std::vector<LocalSldProjectiveHeightSchurEntry> entries;
    std::vector<LocalSldProjectiveHeightSchurGapRow> gaps;
    std::vector<LocalSldProjectiveHeightComponentSchurRow> components;
    SpectralReal maximum_symmetric_geometric_ratio = 0.0L;
    SpectralReal maximum_weighted_row_sum = 0.0L;
    SpectralReal total_component_envelope = 0.0L;
    SpectralReal diagonal_component_envelope = 0.0L;
    SpectralReal weighted_schur_upper_bound = 0.0L;
    SpectralReal upper_bound_ratio = 0.0L;
    SpectralReal commutator_paired_total_envelope = 0.0L;
    SpectralReal commutator_paired_diagonal_envelope = 0.0L;
    SpectralReal commutator_paired_maximum_weighted_row_sum = 0.0L;
    SpectralReal commutator_paired_weighted_schur_upper_bound = 0.0L;
    SpectralReal commutator_paired_upper_bound_ratio = 0.0L;
    SpectralReal commutator_paired_outer_weight = 0.0L;
    SpectralReal commutator_paired_outer_maximum_weighted_row_sum = 0.0L;
    SpectralReal commutator_paired_outer_weighted_schur_upper_bound = 0.0L;
    SpectralReal commutator_paired_outer_upper_bound_ratio = 0.0L;
    std::size_t unscaled_off_diagonal_pair_count = 0;
    std::size_t commutator_paired_unscaled_off_diagonal_pair_count = 0;
    std::size_t
        commutator_paired_outer_unscaled_off_diagonal_pair_count = 0;
    bool finite_matrix_exact = false;
    bool finite_schur_inequality_verified = false;
    bool finite_commutator_paired_schur_inequality_verified = false;
    bool finite_commutator_paired_outer_schur_inequality_verified = false;
    bool cutoff_uniform_weighted_schur_bound_proved = false;
};

class LocalSldProjectiveHeightSchurSummary {
public:
    [[nodiscard]] static LocalSldProjectiveHeightSchurReport summarize(
        const LocalSldProjectiveHeightMatrixReport& matrix);
};

}  // namespace lemma
