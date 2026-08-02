#pragma once

#include "local_sld_projective_height_matrix.hpp"

#include <vector>

namespace lemma {

struct LocalSldProjectiveHeightTailRow {
    int last_core_shell = 0;
    SpectralInteger core_maximum_height = 0;
    SpectralReal core_internal_power_one = 0.0L;
    SpectralReal core_tail_power_one = 0.0L;
    SpectralReal tail_internal_power_one = 0.0L;
    SpectralReal open_power_one = 0.0L;
    SpectralReal open_outer_square_power_one = 0.0L;
    SpectralReal open_advected_commutator_power_one = 0.0L;
    SpectralReal open_advecting_nested_power_one = 0.0L;
    SpectralReal open_enstrophy_normalization_power_one = 0.0L;
    SpectralReal open_palinstrophy_normalization_power_one = 0.0L;
    SpectralReal reconstructed_power_one = 0.0L;
    SpectralReal open_absolute_power_one_sum = 0.0L;
    SpectralReal open_effective_height_pairs = 0.0L;
    SpectralReal dominant_open_pair_fraction = 0.0L;
    SpectralReal open_signed_alignment = 0.0L;
    SpectralReal reconstruction_error = 0.0L;
    SpectralReal component_reconstruction_error = 0.0L;
};

struct LocalSldProjectiveHeightTailReport {
    std::vector<LocalSldProjectiveHeightTailRow> rows;
    SpectralReal maximum_reconstruction_error = 0.0L;
    SpectralReal maximum_component_reconstruction_error = 0.0L;
    bool exact_cumulative_decomposition = false;
    bool uniform_weighted_tail_bound_proved = false;
};

class LocalSldProjectiveHeightTailSummary {
public:
    [[nodiscard]] static LocalSldProjectiveHeightTailReport summarize(
        const LocalSldProjectiveHeightMatrixReport& matrix);
};

}  // namespace lemma
