#include "local_quartic_closure_target.hpp"

#include <algorithm>
#include <cmath>

namespace lemma {
namespace {

SpectralReal relative_error(SpectralReal value, SpectralReal reference) {
    const SpectralReal scale = std::max(
        {std::abs(value), std::abs(reference), 1.0e-30L});
    return std::abs(value - reference) / scale;
}

}  // namespace

LocalQuarticClosureTargetReport LocalQuarticClosureTarget::evaluate(
    const ShiftedCriticalDensityDiagnostic& diagnostic,
    const LocalCriticalDerivativeLedgerReport& derivative_ledger,
    const LocalQuarticReducedReport& reduced) {
    LocalQuarticClosureTargetReport result;
    const SpectralReal stretching = derivative_ledger.signed_stretching;
    const SpectralReal enstrophy = derivative_ledger.enstrophy;
    const SpectralReal palinstrophy = derivative_ledger.palinstrophy;
    const SpectralReal shift = diagnostic.initial_ep_shift;
    const SpectralReal initial_frequency = diagnostic.initial_frequency;
    result.signed_two_entry_bracket = reduced.reduced_pairing +
        derivative_ledger.stretching_local_nonlinear_roles
            .advecting_slot;
    result.absolute_two_entry_bracket =
        std::abs(result.signed_two_entry_bracket);
    if (enstrophy > 0.0L && palinstrophy > 0.0L && shift > 0.0L &&
        initial_frequency > 0.0L) {
        result.candidate_scale = initial_frequency *
            std::pow(shift, 0.25L) *
            std::pow(enstrophy, 1.25L) *
            std::pow(palinstrophy, 0.75L);
    }
    if (result.candidate_scale > 0.0L) {
        result.required_constant_ratio =
            result.absolute_two_entry_bracket / result.candidate_scale;
    }
    const SpectralReal stretching2 = stretching * stretching;
    const SpectralReal stretching4 = stretching2 * stretching2;
    const SpectralReal enstrophy2 = enstrophy * enstrophy;
    const SpectralReal enstrophy3 = enstrophy2 * enstrophy;
    const SpectralReal palinstrophy2 = palinstrophy * palinstrophy;
    const SpectralReal palinstrophy4 = palinstrophy2 * palinstrophy2;
    result.sld_first_denominator_term =
        stretching4 * enstrophy2 * palinstrophy;
    result.sld_shift_denominator_term =
        shift * enstrophy3 * palinstrophy4;
    if (result.sld_first_denominator_term >= 0.0L &&
        result.sld_shift_denominator_term >= 0.0L) {
        result.weighted_geometric_mean = std::pow(
            result.sld_first_denominator_term, 0.75L) *
            std::pow(result.sld_shift_denominator_term, 0.25L);
    }
    result.local_absolute_polynomial_numerator =
        4.0L * std::abs(stretching * stretching2) *
        enstrophy * palinstrophy *
        result.absolute_two_entry_bracket;
    result.reconstructed_scaled_geometric_mean =
        4.0L * initial_frequency *
        result.required_constant_ratio *
        result.weighted_geometric_mean;
    result.geometric_identity_error = relative_error(
        result.local_absolute_polynomial_numerator,
        result.reconstructed_scaled_geometric_mean);
    result.young_upper_bound = initial_frequency *
        result.required_constant_ratio *
        (3.0L * result.sld_first_denominator_term +
         result.sld_shift_denominator_term);
    if (result.young_upper_bound > 0.0L) {
        result.young_ratio =
            result.local_absolute_polynomial_numerator /
            result.young_upper_bound;
    }
    result.certified_local_constant =
        3.0L * result.required_constant_ratio;
    result.algebra_certified =
        result.geometric_identity_error < 1.0e-12L &&
        result.young_ratio <= 1.0L + 1.0e-15L;
    result.finite = reduced.finite && result.algebra_certified &&
        std::isfinite(result.required_constant_ratio);
    return result;
}

}  // namespace lemma
