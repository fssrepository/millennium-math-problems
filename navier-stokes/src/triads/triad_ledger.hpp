#pragma once

#include "triad_partition.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace lemma {

enum class TriadLowRole : std::size_t {
    advecting = 0,
    advected = 1,
    target = 2,
    tied = 3
};

struct TriadGapLedgerRow {
    int dyadic_gap = 0;
    std::size_t interactions = 0;
    SpectralReal signed_stretching = 0.0L;
    SpectralReal absolute_pair_stretching = 0.0L;
    std::array<std::size_t, 4> interactions_by_low_role{};
    std::array<SpectralReal, 4> signed_stretching_by_low_role{};
    std::array<SpectralReal, 4> absolute_stretching_by_low_role{};
};

struct TriadLedgerReport {
    SpectralReal signed_total = 0.0L;
    SpectralReal signed_local = 0.0L;
    SpectralReal signed_nonlocal = 0.0L;
    SpectralReal absolute_pair_total = 0.0L;
    SpectralReal partition_residual = 0.0L;
    std::vector<TriadGapLedgerRow> gaps;
};

class TriadLedger {
public:
    [[nodiscard]] static TriadLedgerReport analyze(
        const SpectralState& state);
    [[nodiscard]] static const char* low_role_name(TriadLowRole role);
};

}  // namespace lemma
