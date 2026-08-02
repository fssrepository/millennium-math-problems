#pragma once

#include "helical_triad_ledger.hpp"

#include <cstdint>

namespace lemma {

enum class HelicalLocalSpread { all, equal, narrow, broad };

struct HelicalSectorSelection {
    std::uint8_t sector_mask = UINT8_C(0xff);
    HelicalLocalSpread spread = HelicalLocalSpread::all;

    [[nodiscard]] static constexpr HelicalSectorSelection homochiral() {
        return {static_cast<std::uint8_t>((UINT8_C(1) << 0U) |
                                          (UINT8_C(1) << 7U))};
    }
    [[nodiscard]] static constexpr HelicalSectorSelection heterochiral() {
        return {static_cast<std::uint8_t>(
            UINT8_C(0xff) ^ homochiral().sector_mask)};
    }
    [[nodiscard]] constexpr bool includes(std::size_t sector) const {
        return (sector_mask & static_cast<std::uint8_t>(
                                  UINT8_C(1) << sector)) != 0U;
    }
    [[nodiscard]] constexpr HelicalSectorSelection with_spread(
        HelicalLocalSpread value) const {
        return {sector_mask, value};
    }
    [[nodiscard]] bool includes_spread(
        WaveVector first, WaveVector second, WaveVector third) const;
};

struct HelicalSectorObjectiveValue {
    SpectralReal signed_local_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal critical_integrand = 0.0L;
};

class HelicalSectorObjective {
public:
    [[nodiscard]] static HelicalSectorObjectiveValue evaluate(
        const SpectralState& state, HelicalSectorSelection selection);
    [[nodiscard]] static SpectralIncrement signed_stretching_gradient(
        const SpectralState& state, HelicalSectorSelection selection);
    [[nodiscard]] static SpectralIncrement critical_integrand_gradient(
        const SpectralState& state, HelicalSectorSelection selection);
};

}  // namespace lemma
