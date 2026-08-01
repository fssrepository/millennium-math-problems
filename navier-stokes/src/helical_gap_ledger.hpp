#pragma once

#include "helical_triad_ledger.hpp"

#include <cstddef>
#include <vector>

namespace lemma {

struct HelicalGapLedgerRow {
    int dyadic_gap = 0;
    std::size_t interactions = 0;
    SpectralReal homochiral_signed = 0.0L;
    SpectralReal homochiral_absolute = 0.0L;
    SpectralReal heterochiral_signed = 0.0L;
    SpectralReal heterochiral_absolute = 0.0L;
};

struct HelicalGapLedgerReport {
    SpectralReal signed_total = 0.0L;
    SpectralReal homochiral_signed_total = 0.0L;
    SpectralReal heterochiral_signed_total = 0.0L;
    SpectralReal maximum_gap_reconstruction_residual = 0.0L;
    SpectralReal total_reconstruction_residual = 0.0L;
    std::vector<HelicalGapLedgerRow> gaps;
};

class HelicalGapLedger {
public:
    [[nodiscard]] static HelicalGapLedgerReport analyze(
        const SpectralState& state);
};

}  // namespace lemma
