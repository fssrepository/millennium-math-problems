#include "transition_block_scaling.hpp"

#include <cmath>
#include <stdexcept>

namespace lemma {

TransitionBlockScalingReport TransitionBlockScaling::analyze() {
    TransitionBlockScalingReport report;
    report.post_young_logarithm_power =
        report.young_remainder_conjugate *
        report.band_count_logarithm_power;
    report.post_young_enstrophy_power =
        report.young_remainder_conjugate *
        report.vortex_enstrophy_power;
    report.required_pointwise_depletion_power =
        (report.post_young_enstrophy_power -
         report.energy_time_integrable_enstrophy_power) /
        report.young_remainder_conjugate;
    report.logarithmic_band_count_changes_polynomial_power = false;
    report.energy_identity_closes_transition_block =
        report.post_young_enstrophy_power <=
            report.energy_time_integrable_enstrophy_power;
    return report;
}

long double TransitionBlockScaling::normalized_remainder(long double z) {
    if (!(z >= 0.0L) || !std::isfinite(z)) {
        throw std::invalid_argument("invalid transition-block enstrophy");
    }
    const long double logarithm = std::log1p(z);
    return z * z * logarithm * logarithm * logarithm * logarithm;
}

}  // namespace lemma
