#pragma once

#include <cstdint>
#include <vector>

namespace lemma {

using ShellGeometryReal = long double;

struct PeriodicShellCountRow {
    int shell = 0;
    std::uint64_t exact_lattice_modes = 0;
    std::uint64_t cube_upper_bound = 0;
    ShellGeometryReal normalized_count_ratio = 0.0L;
};

struct PeriodicShellGeometryCertificate {
    int maximum_enumerated_shell = 0;
    int overlap_samples = 0;
    std::uint64_t seed = 0;
    std::vector<PeriodicShellCountRow> shell_counts;
    ShellGeometryReal lattice_count_constant = 64.0L;
    ShellGeometryReal l2_to_linf_bernstein_constant = 8.0L;
    ShellGeometryReal gradient_bernstein_constant = 16.0L;
    int separated_high_shell_neighbor_width = 1;
    ShellGeometryReal one_gain_overlap_constant = 0.0L;
    ShellGeometryReal three_gain_overlap_constant = 0.0L;
    ShellGeometryReal low_advecting_frequency_constant = 16.0L;
    ShellGeometryReal low_advected_frequency_constant = 8.0L;
    ShellGeometryReal low_target_frequency_constant = 8.0L;
    ShellGeometryReal low_advecting_block_constant = 128.0L;
    ShellGeometryReal low_advected_block_constant = 64.0L;
    ShellGeometryReal low_target_block_constant = 64.0L;
    ShellGeometryReal ft1_one_gain_constant = 0.0L;
    ShellGeometryReal ft1_three_gain_constant = 0.0L;
    ShellGeometryReal maximum_count_ratio = 0.0L;
    ShellGeometryReal maximum_one_gain_overlap_ratio = 0.0L;
    ShellGeometryReal maximum_three_gain_overlap_ratio = 0.0L;
    bool separated_high_shells_are_adjacent = true;
    bool all_bounds_hold = true;
};

class PeriodicShellGeometry {
public:
    // Hard shells are S_j={k in Z^3: 2^j<=|k|<2^(j+1)}, j>=0,
    // on the 2*pi torus with normalized Haar measure.
    [[nodiscard]] static PeriodicShellGeometryCertificate certify(
        int maximum_enumerated_shell, int overlap_samples,
        std::uint64_t seed);
};

}  // namespace lemma
