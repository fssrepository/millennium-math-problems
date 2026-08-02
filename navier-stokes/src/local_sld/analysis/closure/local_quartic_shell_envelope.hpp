#pragma once

#include "local_quartic_shell_ledger.hpp"

#include <vector>

namespace lemma {

struct LocalQuarticShellEnvelopeRow {
    int shell = 0;
    SpectralReal radius = 0.0L;
    SpectralReal neighborhood_energy = 0.0L;
    SpectralReal actual_local_advection_h1_squared = 0.0L;
    SpectralReal explicit_bound = 0.0L;
    SpectralReal bound_ratio = 0.0L;
};

struct LocalQuarticShellEnvelopeReport {
    std::vector<LocalQuarticShellEnvelopeRow> shells;
    SpectralReal explicit_constant = 46656.0L;
    int frequency_power = 7;
    int required_envelope_power = 3;
    int palinstrophy_power = 4;
    int residual_frequency_gain = -1;
    SpectralReal maximum_bound_ratio = 0.0L;
    SpectralReal maximum_target_frequency_ratio = 0.0L;
    SpectralReal maximum_advected_frequency_ratio = 0.0L;
    SpectralReal maximum_interaction_count_ratio = 0.0L;
    SpectralReal neighborhood_h3_moment = 0.0L;
    SpectralReal neighborhood_h4_moment = 0.0L;
    SpectralReal shell_product_sum = 0.0L;
    SpectralReal state_enstrophy = 0.0L;
    SpectralReal state_palinstrophy = 0.0L;
    SpectralReal interpolated_h3_bound = 0.0L;
    SpectralReal h3_overlap_constant = 73.0L / 8.0L;
    SpectralReal h4_overlap_constant = 273.0L / 16.0L;
    SpectralReal actual_global_local_advection_h1_squared = 0.0L;
    SpectralReal global_shell_product_bound = 0.0L;
    SpectralReal global_zp_bound = 0.0L;
    SpectralReal global_shell_product_ratio = 0.0L;
    SpectralReal global_zp_ratio = 0.0L;
    bool cutoff_independent = true;
    bool all_inputs_in_neighboring_shells = false;
    bool all_geometry_checks_hold = false;
    bool global_summation_holds = false;
    bool all_bounds_hold = false;
};

class LocalQuarticShellEnvelope {
public:
    // For R=2^j and energy in j-1..j+1, certifies
    // ||A^(1/2) B_local||_{2,j}^2 <= 46656 R^7 E_near^2.
    // Summing R^7 E_near^2=(R^3 E_near)(R^4 E_near) and using
    // H^(3/2)^2 <= sqrt(Z P) gives a cutoff-independent global bound
    // ||A^(1/2)B_local||_2^2 <= C sqrt(Z) P^(3/2).
    [[nodiscard]] static LocalQuarticShellEnvelopeReport analyze(
        const SpectralState& state,
        const LocalQuarticShellReport& shell_ledger);
};

}  // namespace lemma
