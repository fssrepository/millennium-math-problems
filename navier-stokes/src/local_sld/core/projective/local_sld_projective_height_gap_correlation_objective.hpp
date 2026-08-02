#pragma once

#include "projective_advection_decomposition.hpp"
#include "triad_partition.hpp"

namespace lemma {

struct LocalSldProjectiveHeightGapCorrelationObjectiveValue {
    int first_shell = 0;
    int second_shell = 0;
    int shell_gap = 0;
    SpectralReal gram_pairing = 0.0L;
    SpectralReal first_h2_norm2 = 0.0L;
    SpectralReal second_h2_norm2 = 0.0L;
    SpectralReal correlation_squared = 0.0L;
    SpectralReal half_decay_weighted_correlation = 0.0L;
    SpectralReal weighted_correlation_squared = 0.0L;
    bool finite = false;
};

// Exact-gradient adversary objective for the proposed PNT-13 height-gap
// almost-orthogonality estimate. The optimized value is
// 4^gap |<A b_i,A b_j>|^2/(||A b_i||^2 ||A b_j||^2).
class LocalSldProjectiveHeightGapCorrelationObjective {
public:
    LocalSldProjectiveHeightGapCorrelationObjective(
        TriadSelection selection,
        int first_shell,
        int second_shell,
        int threads = 12);

    [[nodiscard]]
    LocalSldProjectiveHeightGapCorrelationObjectiveValue evaluate(
        const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    TriadSelection selection_;
    int first_shell_ = 0;
    int second_shell_ = 0;
    int threads_ = 12;
};

}  // namespace lemma
