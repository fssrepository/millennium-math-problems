#pragma once

#include "proof_scaling.hpp"
#include "spectral_state.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct EqualLowTriadCountRow {
    int cutoff = 0;
    std::size_t ordered_pairs = 0;
    std::size_t maximum_input_degree = 0;
    std::size_t maximum_target_degree = 0;
    std::size_t elementary_degree_bound = 0;
};

struct EqualLowTriadGeometryCertificate {
    int target_length_multiplier = 3;
    int maximum_cutoff = 0;
    SpectralReal maximum_input_degree_ratio = 0.0L;
    SpectralReal maximum_target_degree_ratio = 0.0L;
    bool exact_fixed_angle_relation = false;
    bool all_degree_bounds_hold = true;
    std::vector<EqualLowTriadCountRow> rows;
};

class EqualLowTriadGeometry {
public:
    [[nodiscard]] static EqualLowTriadGeometryCertificate certify(
        int maximum_cutoff,
        int target_length_multiplier);
};

}  // namespace lemma
