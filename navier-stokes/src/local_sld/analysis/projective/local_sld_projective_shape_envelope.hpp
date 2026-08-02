#pragma once

#include "local_sld_remainder_projective_ledger.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct LocalSldProjectiveShapeShellRow {
    int dyadic_height_level = 0;
    SpectralInteger minimum_primitive_height = 0;
    SpectralInteger maximum_primitive_height = 0;
    std::size_t shape_count = 0;
    SpectralReal signed_power_one_total = 0.0L;
    SpectralReal absolute_power_one_sum = 0.0L;
    SpectralReal squared_power_one_sum = 0.0L;
    SpectralReal effective_shapes = 0.0L;
    SpectralReal signed_alignment = 0.0L;
    SpectralReal absolute_fraction = 0.0L;
    SpectralReal weighted_angle_sine_squared = 0.0L;
    SpectralReal weighted_length_aspect_ratio = 0.0L;
};

struct LocalSldProjectiveShapeEnvelopeReport {
    std::size_t projective_shape_count = 0;
    SpectralReal expected_signed_total = 0.0L;
    SpectralReal reconstructed_signed_total = 0.0L;
    SpectralReal reconstruction_error = 0.0L;
    SpectralReal absolute_total = 0.0L;
    SpectralReal contribution_weighted_angle_sine_squared = 0.0L;
    SpectralReal contribution_weighted_length_aspect_ratio = 0.0L;
    SpectralReal contribution_weighted_dyadic_span = 0.0L;
    SpectralReal primitive_height_half_moment = 0.0L;
    SpectralReal primitive_height_first_moment = 0.0L;
    SpectralReal fitted_absolute_height_shell_slope = 0.0L;
    bool exact_reconstruction = false;
    bool summable_projective_shape_envelope_proved = false;
    std::vector<LocalSldProjectiveShapeShellRow> height_shells;
};

class LocalSldProjectiveShapeEnvelope {
public:
    [[nodiscard]] static LocalSldProjectiveShapeEnvelopeReport analyze(
        const LocalSldRemainderProjectiveReport& projective);
};

}  // namespace lemma
