#pragma once

#include "proof_scaling.hpp"
#include "spectral_state.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace lemma {

struct ProjectiveTriadCountRow {
    int cutoff = 0;
    std::size_t ordered_pairs = 0;
    std::size_t maximum_input_degree = 0;
    std::size_t maximum_target_degree = 0;
    std::size_t elementary_degree_bound = 0;
};

struct ProjectiveTriadGeometryCertificate {
    std::array<SpectralInteger, 3> primitive_squared_lengths{};
    int maximum_cutoff = 0;
    std::size_t ordered_role_count = 0;
    SpectralReal maximum_input_degree_ratio = 0.0L;
    SpectralReal maximum_target_degree_ratio = 0.0L;
    bool primitive_signature = false;
    bool triangle_feasible = false;
    bool fixed_plane_sphere_geometry = false;
    bool all_degree_bounds_hold = true;
    std::vector<ProjectiveTriadCountRow> rows;
};

class ProjectiveTriadGeometry {
public:
    [[nodiscard]] static ProjectiveTriadGeometryCertificate certify(
        int maximum_cutoff,
        std::array<SpectralInteger, 3> primitive_squared_lengths);
};

}  // namespace lemma
