#include "local_signature_density.hpp"

#include "local_signature_objective.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace lemma {

LocalSignatureDensitySample LocalSignatureDensity::evaluate(
    const SpectralState& state) {
    const LocalSignatureObjectiveValue signature =
        LocalSignatureObjective::evaluate(state);
    LocalSignatureDensitySample result;
    result.amplification = signature.signed_amplification;
    result.amplification_fourth = std::pow(
        result.amplification, 4.0L);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal energy = std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
        result.enstrophy += wave2 * energy;
        result.palinstrophy += wave2 * wave2 * energy;
    }
    if (!(result.enstrophy > 0.0L) ||
        !(result.palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal denominator = result.enstrophy *
        result.palinstrophy * result.palinstrophy * result.palinstrophy;
    const SpectralReal transfer2 = signature.signed_local_transfer *
        signature.signed_local_transfer;
    result.critical_density = transfer2 * transfer2 / denominator;
    result.square_signature_density =
        signature.squared_signature_transfer *
        signature.squared_signature_transfer / denominator;
    const SpectralReal factored = result.amplification_fourth *
        result.square_signature_density;
    result.factorization_residual = std::abs(
        result.critical_density - factored) /
        std::max({std::numeric_limits<SpectralReal>::min(),
                  std::abs(result.critical_density), std::abs(factored)});
    return result;
}

}  // namespace lemma
