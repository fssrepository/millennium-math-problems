#pragma once

#include "spectral_state.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct ProjectiveFanGeometryRow {
    int cutoff = 0;
    std::size_t pair_count = 0;
    std::size_t primitive_shape_count = 0;
    SpectralReal target_synthesis_ratio = 0.0L;
    SpectralReal synthesis_ratio_per_pair = 0.0L;
    SpectralReal quadratic_lower_bound_ratio = 0.0L;
    bool every_shape_unique = false;
    bool quadratic_lower_bound_verified = false;
};

struct ProjectiveFanGeometryCertificate {
    int maximum_cutoff = 0;
    bool local_frequency_geometry_proved = false;
    bool target_is_unique_largest_role_proved = false;
    bool primitive_shape_injectivity_proved = false;
    bool aligned_target_coefficients_proved = false;
    bool pair_count_quadratic_lower_bound_proved = false;
    bool target_synthesis_quadratic_lower_bound_proved = false;
    bool target_synthesis_unbounded_proved = false;
    bool stretching_support_disjointness_proved = false;
    bool exact_zero_stretching_proved = false;
    bool exact_zero_power_one_product_proved = false;
    std::vector<ProjectiveFanGeometryRow> rows;
};

class ProjectiveFanGeometry {
public:
    [[nodiscard]] static ProjectiveFanGeometryCertificate certify(
        int maximum_cutoff);
};

}  // namespace lemma
