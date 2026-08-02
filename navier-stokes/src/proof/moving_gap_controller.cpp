#include "moving_gap_controller.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {

MovingGapDecision MovingGapController::decide(
    SpectralReal enstrophy, int base_gap) {
    if (!(enstrophy >= 0.0L) || !std::isfinite(enstrophy) ||
        base_gap < 0 || base_gap > 1024) {
        throw std::invalid_argument("invalid moving-gap parameters");
    }
    MovingGapDecision result;
    result.base_gap = base_gap;
    result.enstrophy = enstrophy;
    if (enstrophy > 1.0L) {
        result.logarithmic_gap = static_cast<int>(
            std::ceil(std::log2(enstrophy)));
    }
    result.minimum_gap = result.base_gap + result.logarithmic_gap;
    if (enstrophy > 0.0L) {
        result.normalized_cubic_remainder_ratio = std::exp2(
            2.0L * std::log2(enstrophy) -
            2.0L * static_cast<SpectralReal>(result.logarithmic_gap));
        result.base_weighted_remainder_ratio = std::scalbn(
            result.normalized_cubic_remainder_ratio,
            -2 * result.base_gap);
    }
    return result;
}

}  // namespace lemma
