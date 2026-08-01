#pragma once

#include "triad_partition.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct TriadCommutatorGapRow {
    int dyadic_gap = 0;
    std::size_t pairs = 0;
    SpectralReal signed_paired_stretching = 0.0L;
    SpectralReal absolute_unpaired_stretching = 0.0L;
    SpectralReal absolute_paired_stretching = 0.0L;
    SpectralReal relative_unweighted_cancellation_residual = 0.0L;
    SpectralReal relative_weighted_identity_residual = 0.0L;
    SpectralReal maximum_frequency_gain_ratio = 0.0L;
};

struct TriadCommutatorReport {
    std::size_t pairs = 0;
    SpectralReal signed_paired_stretching = 0.0L;
    SpectralReal absolute_unpaired_stretching = 0.0L;
    SpectralReal absolute_paired_stretching = 0.0L;
    SpectralReal relative_unweighted_cancellation_residual = 0.0L;
    SpectralReal relative_weighted_identity_residual = 0.0L;
    SpectralReal maximum_frequency_gain_ratio = 0.0L;
    std::vector<TriadCommutatorGapRow> gaps;
};

class TriadCommutator {
public:
    // Pairs (p,q,k) with (p,-k,-q) when p is the unique low wave.
    // Incompressibility and reality give T(p,q,k)+T(p,-k,-q)=0
    // before enstrophy weighting, leaving (|k|^2-|q|^2)T(p,q,k).
    [[nodiscard]] static TriadCommutatorReport analyze(
        const SpectralState& state);
};

}  // namespace lemma
