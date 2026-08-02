#pragma once

#include <cstdint>
#include <vector>

namespace lemma {

using ShellReal = long double;

struct DyadicShellBoundReport {
    int shells = 0;
    int minimum_gap = 0;
    ShellReal enstrophy = 0.0L;
    ShellReal palinstrophy = 0.0L;
    ShellReal high_moment = 0.0L;
    ShellReal high_moment_bound = 0.0L;
    ShellReal high_moment_ratio = 0.0L;
    ShellReal low_one_derivative_constant = 0.0L;
    ShellReal low_three_derivative_constant = 0.0L;
    ShellReal maximum_low_one_derivative_ratio = 0.0L;
    ShellReal maximum_low_three_derivative_ratio = 0.0L;
    ShellReal one_gain_tail = 0.0L;
    ShellReal one_gain_tail_bound = 0.0L;
    ShellReal one_gain_tail_ratio = 0.0L;
    ShellReal three_gain_tail = 0.0L;
    ShellReal three_gain_tail_bound = 0.0L;
    ShellReal three_gain_tail_ratio = 0.0L;
};

struct DyadicShellRandomCertificate {
    int shells = 0;
    int minimum_gap = 0;
    int samples = 0;
    std::uint64_t seed = 0;
    ShellReal maximum_high_moment_ratio = 0.0L;
    ShellReal maximum_low_one_derivative_ratio = 0.0L;
    ShellReal maximum_low_three_derivative_ratio = 0.0L;
    ShellReal maximum_one_gain_tail_ratio = 0.0L;
    ShellReal maximum_three_gain_tail_ratio = 0.0L;
    bool all_bounds_hold = true;
};

class DyadicShellBounds {
public:
    // a_j is the nonnegative L2 amplitude of shell j, j>=0. The model norms
    // are Z=sum 2^(2j)a_j^2 and P=sum 2^(4j)a_j^2.
    [[nodiscard]] static DyadicShellBoundReport analyze(
        const std::vector<ShellReal>& amplitudes, int minimum_gap);

    // Deterministic numerical regression for the cutoff-independent sequence
    // inequalities. It does not certify Bernstein or LP projection constants.
    [[nodiscard]] static DyadicShellRandomCertificate verify_random(
        int shells, int minimum_gap, int samples, std::uint64_t seed);
};

}  // namespace lemma
