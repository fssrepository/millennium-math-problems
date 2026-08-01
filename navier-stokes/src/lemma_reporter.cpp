#include "lemma_reporter.hpp"

#include <iomanip>
#include <ostream>

namespace lemma {

void LemmaReporter::write_console(const LemmaReport& report, std::ostream& out) {
    out << std::setprecision(10)
        << "L2 scaling certificate\n"
        << "  admissible absorbable exponent samples: " << report.candidate_count << '\n'
        << "  minimum post-Young Z power:             "
        << report.minimum_young_power << '\n'
        << "  minimizer (a,b,c):                      ("
        << report.minimizer_energy << ',' << report.minimizer_enstrophy << ','
        << report.minimizer_palinstrophy << ")\n"
        << "  Young multiplier power r:               "
        << report.young_multiplier_power << '\n'
        << "  pointwise-linear depletion delta:        "
        << report.pointwise_depletion_power << '\n'
        << "  energy-integrable depletion delta:       "
        << report.energy_depletion_power << '\n'
        << "  universal quarter-depletion identity:    "
        << (report.universal_quarter_depletion ? "verified" : "FAILED") << '\n'
        << "  globally closing E-Z-P candidate:       "
        << (report.closing_candidate_exists ? "YES" : "none") << "\n\n"
        << "Concentration scaling\n"
        << "  fixed-energy pointwise Q exponent:       "
        << report.fixed_energy_q_exponent << " (must be <= 0)\n"
        << "  pointwise Q candidate:                   "
        << (report.pointwise_q_scale_compatible ? "survives" : "REJECTED") << '\n'
        << "  integrated D^4 Z^2 dt exponent:          "
        << report.integrated_l4_exponent << '\n'
        << "  time-integrated L4-A:                    "
        << (report.integrated_l4_scale_critical ? "scale-critical" : "REJECTED")
        << "\n  strong-Q factorization D^4 Z^2=Q Z:     "
        << (report.exact_strong_l4_factorization ? "verified" : "FAILED")
        << "\n  uniform trajectory Q closes L4 via E:    "
        << (report.uniform_q_closes_l4 ? "verified" : "FAILED") << "\n\n"
        << "Dyadic far-tail scaling\n"
        << "  advecting/advected gap decay:           2^(-"
        << report.dyadic_advecting_gap_decay << " m)\n"
        << "  target gap decay:                       2^(-"
        << report.dyadic_target_gap_decay << " m)\n"
        << "  fourth-power density gap decay:         2^(-"
        << report.dyadic_l4_density_gap_decay << " m)\n"
        << "  remaining density Z power:              "
        << report.dyadic_l4_density_enstrophy_power << '\n'
        << "  frequency-gap series:                   "
        << (report.dyadic_tail_summable ? "summable" : "NOT SUMMABLE")
        << "\n  closes in time from energy identity:     "
        << (report.dyadic_energy_closes_time_integral ? "YES" : "no")
        << "\n  post-Young remainder:                    2^(-"
        << report.dyadic_post_young_gap_decay << " m) Z^"
        << report.dyadic_post_young_enstrophy_power
        << "\n  moving gap m~log2(Z) leaves Z power:     "
        << report.moving_gap_remaining_enstrophy_power
        << "\n  moving gap closes geometric far tail:    "
        << (report.moving_gap_closes_far_tail ? "verified" : "FAILED")
        << "\n  shell-sequence samples:                  "
        << report.dyadic_shell_bounds.samples << " x "
        << report.dyadic_shell_bounds.shells
        << " shells\n  max high-moment Holder ratio:           "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_high_moment_ratio)
        << "\n  max one-gain tail ratio:                 "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_one_gain_tail_ratio)
        << "\n  max three-gain tail ratio:               "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_three_gain_tail_ratio)
        << "\n  shell-sequence bounds:                   "
        << (report.dyadic_shell_bounds.all_bounds_hold
                ? "verified"
                : "FAILED")
        << "\n  hard-shell lattice count constant:       "
        << static_cast<double>(
               report.periodic_shell_geometry.lattice_count_constant)
        << "\n  hard-shell L2-to-Linf constant:          "
        << static_cast<double>(
               report.periodic_shell_geometry
                   .l2_to_linf_bernstein_constant)
        << "\n  high-shell 5/2-overlap constants:        "
        << static_cast<double>(
               report.periodic_shell_geometry.one_gain_overlap_constant)
        << " / "
        << static_cast<double>(
               report.periodic_shell_geometry.three_gain_overlap_constant)
        << "\n  explicit FT-1 constants (one/three):     "
        << static_cast<double>(
               report.periodic_shell_geometry.ft1_one_gain_constant)
        << " / "
        << static_cast<double>(
               report.periodic_shell_geometry.ft1_three_gain_constant)
        << "\n  periodic shell geometry:                "
        << (report.periodic_shell_geometry.all_bounds_hold
                ? "verified"
                : "FAILED")
        << "\n  direct FT-1 Fourier samples / max ratio: "
        << report.periodic_tail_bound.samples << " / "
        << static_cast<double>(
               report.periodic_tail_bound.maximum_bound_ratio)
        << "\n  explicit periodic FT-1 regression:       "
        << (report.periodic_tail_bound.all_bounds_hold
                ? "verified"
                : "FAILED")
        << "\n  moving-tail Young samples / max ratio:   "
        << report.far_tail_closure.samples << " / "
        << static_cast<double>(
               report.far_tail_closure.maximum_normalized_remainder_ratio)
        << "\n  cutoff-independent moving-tail closure: "
        << (report.far_tail_closure.all_bounds_hold
                ? "verified"
                : "FAILED")
        << "\n  transition remainder after Young:        log(1+Z)^"
        << report.transition_block_scaling.post_young_logarithm_power.str()
        << " Z^"
        << report.transition_block_scaling.post_young_enstrophy_power.str()
        << "\n  required local pointwise depletion:      Z^(-"
        << report.transition_block_scaling
               .required_pointwise_depletion_power.str()
        << ")\n  log-band count closes transition block:  "
        << (report.transition_block_scaling
                    .energy_identity_closes_transition_block
                ? "YES"
                : "no")
        << "\n  helical sector reconstruction residual: "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_local_reconstruction_residual)
        << "\n  max pure homochiral local stretching:    "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_pure_homochiral_local_stretching)
        << "\n  homochiral zero-stretching candidate:    "
        << (report.helical_triad_certificate
                    .nonzero_pure_homochiral_local_seen
                ? "REJECTED"
                : "not rejected")
        << "\n\n"
        << "Fourier-Galerkin triad checks\n"
        << "  modes / samples:                        " << report.triad_modes
        << " / " << report.triad_samples << '\n'
        << "  max normalized energy residual:         "
        << static_cast<double>(report.energy_residual) << '\n'
        << "  max divergence residual:                "
        << static_cast<double>(report.divergence_residual) << '\n'
        << "  max detailed-triad residual:            "
        << static_cast<double>(report.detailed_triad_residual) << '\n'
        << "  nonzero vortex stretching:              "
        << (report.nonzero_vortex_stretching ? "observed" : "NOT OBSERVED") << '\n'
        << "  max classical scale-invariant ratio:    "
        << static_cast<double>(report.classical_ratio) << '\n'
        << "  max nonlocal absolute fraction:         "
        << static_cast<double>(report.nonlocal_absolute_fraction) << '\n'
        << "  max net-flux interaction efficiency:    "
        << static_cast<double>(report.flux_efficiency) << "\n\n"
        << "Conclusion: the E-Z-P monomial route cannot close the global enstrophy "
           "estimate. The next lemma must exploit an additional dynamic/geometric "
           "quantity while retaining cutoff-independent constants.\n";
}

void LemmaReporter::write_json(const LemmaReport& report, std::ostream& out) {
    out << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-lemma-certificate-v2\",\n"
        << "  \"target\": \"L2 and strong-L4 reduction\",\n"
        << "  \"claim_scope\": \"nonnegative monomial E^a Z^b P^c bounds with NS scaling\",\n"
        << "  \"scaling\": {\n"
        << "    \"amplitude_equation\": \"2(a+b+c)=3\",\n"
        << "    \"space_equation\": \"-a+b+3c=3\",\n"
        << "    \"minimum_post_young_Z_power\": \"" << report.minimum_young_power
        << "\",\n    \"minimizer\": {\"a\": \"" << report.minimizer_energy
        << "\", \"b\": \"" << report.minimizer_enstrophy
        << "\", \"c\": \"" << report.minimizer_palinstrophy
        << "\", \"young_multiplier_power\": \"" << report.young_multiplier_power
        << "\", \"pointwise_linear_depletion_power\": \""
        << report.pointwise_depletion_power
        << "\", \"energy_integrable_depletion_power\": \""
        << report.energy_depletion_power << "\"},\n"
        << "    \"all_absorbable_candidates_require_energy_depletion_power\": \"1/4\",\n"
        << "    \"quarter_depletion_identity_verified\": "
        << (report.universal_quarter_depletion ? "true" : "false") << ",\n"
        << "    \"globally_closing_candidate_found\": "
        << (report.closing_candidate_exists ? "true" : "false") << "\n  },\n"
        << "  \"concentration_scaling\": {\n"
        << "    \"map\": \"u_lambda(x,t)=lambda*u(lambda*x,lambda^2*t)\",\n"
        << "    \"fixed_energy_Q_exponent\": \"" << report.fixed_energy_q_exponent
        << "\",\n    \"pointwise_Q_candidate_scale_compatible\": "
        << (report.pointwise_q_scale_compatible ? "true" : "false") << ",\n"
        << "    \"D4Z2_density_exponent\": \"" << report.critical_density_exponent
        << "\",\n    \"dt_exponent\": \"" << report.time_exponent
        << "\",\n    \"integrated_L4_exponent\": \"" << report.integrated_l4_exponent
        << "\",\n    \"integrated_L4_scale_critical\": "
        << (report.integrated_l4_scale_critical ? "true" : "false") << "\n  },\n"
        << "  \"strong_L4_reduction\": {\n"
        << "    \"Q\": \"D^4 Z\",\n"
        << "    \"critical_density_factorization\": \"D^4 Z^2 = Q Z\",\n"
        << "    \"exact_factorization_verified\": "
        << (report.exact_strong_l4_factorization ? "true" : "false") << ",\n"
        << "    \"energy_identity\": \"2 nu integral Z dt = E(0)-E(T)\",\n"
        << "    \"conditional_bound\": \"integral D^4 Z^2 dt <= (sup Q) E(0)/(2 nu)\",\n"
        << "    \"uniform_Q_closes_L4\": "
        << (report.uniform_q_closes_l4 ? "true" : "false") << "\n  },\n"
        << "  \"dyadic_tail_scaling\": {\n"
        << "    \"low_advecting_gap_decay\": \""
        << report.dyadic_advecting_gap_decay
        << "\",\n    \"low_advected_gap_decay\": \""
        << report.dyadic_advecting_gap_decay
        << "\",\n    \"low_target_gap_decay\": \""
        << report.dyadic_target_gap_decay
        << "\",\n    \"L4_density_gap_decay\": \""
        << report.dyadic_l4_density_gap_decay
        << "\",\n    \"L4_density_enstrophy_power\": \""
        << report.dyadic_l4_density_enstrophy_power
        << "\",\n    \"frequency_tail_summable\": "
        << (report.dyadic_tail_summable ? "true" : "false")
        << ",\n    \"energy_identity_closes_time_integral\": "
        << (report.dyadic_energy_closes_time_integral ? "true" : "false")
        << ",\n    \"post_young_gap_decay\": \""
        << report.dyadic_post_young_gap_decay
        << "\",\n    \"post_young_enstrophy_power\": \""
        << report.dyadic_post_young_enstrophy_power
        << "\",\n    \"moving_gap_log_enstrophy_slope\": \""
        << report.moving_gap_log_enstrophy_slope
        << "\",\n    \"moving_gap_remaining_enstrophy_power\": \""
        << report.moving_gap_remaining_enstrophy_power
        << "\",\n    \"moving_gap_closes_far_tail\": "
        << (report.moving_gap_closes_far_tail ? "true" : "false")
        << "\n  },\n"
        << "  \"dyadic_shell_bounds\": {\n"
        << "    \"shells\": " << report.dyadic_shell_bounds.shells
        << ",\n    \"minimum_gap\": "
        << report.dyadic_shell_bounds.minimum_gap
        << ",\n    \"samples\": " << report.dyadic_shell_bounds.samples
        << ",\n    \"seed\": " << report.dyadic_shell_bounds.seed
        << ",\n    \"low_one_derivative_constant\": \"sqrt(2)\",\n"
        << "    \"low_three_derivative_constant\": \"sqrt(8/7)\",\n"
        << "    \"high_moment_constant\": \"1\",\n"
        << "    \"maximum_high_moment_ratio\": "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_high_moment_ratio)
        << ",\n    \"maximum_low_one_derivative_ratio\": "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_low_one_derivative_ratio)
        << ",\n    \"maximum_low_three_derivative_ratio\": "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_low_three_derivative_ratio)
        << ",\n    \"maximum_one_gain_tail_ratio\": "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_one_gain_tail_ratio)
        << ",\n    \"maximum_three_gain_tail_ratio\": "
        << static_cast<double>(
               report.dyadic_shell_bounds.maximum_three_gain_tail_ratio)
        << ",\n    \"all_bounds_hold\": "
        << (report.dyadic_shell_bounds.all_bounds_hold ? "true" : "false")
        << "\n  },\n"
        << "  \"periodic_shell_geometry\": {\n"
        << "    \"maximum_enumerated_shell\": "
        << report.periodic_shell_geometry.maximum_enumerated_shell
        << ",\n    \"overlap_samples\": "
        << report.periodic_shell_geometry.overlap_samples
        << ",\n    \"seed\": " << report.periodic_shell_geometry.seed
        << ",\n    \"lattice_count_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry.lattice_count_constant)
        << ",\n    \"l2_to_linf_bernstein_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry
                   .l2_to_linf_bernstein_constant)
        << ",\n    \"gradient_bernstein_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry.gradient_bernstein_constant)
        << ",\n    \"separated_high_shell_neighbor_width\": "
        << report.periodic_shell_geometry.separated_high_shell_neighbor_width
        << ",\n    \"one_gain_overlap_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry.one_gain_overlap_constant)
        << ",\n    \"three_gain_overlap_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry.three_gain_overlap_constant)
        << ",\n    \"ft1_one_gain_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry.ft1_one_gain_constant)
        << ",\n    \"ft1_three_gain_constant\": "
        << static_cast<double>(
               report.periodic_shell_geometry.ft1_three_gain_constant)
        << ",\n    \"maximum_count_ratio\": "
        << static_cast<double>(
               report.periodic_shell_geometry.maximum_count_ratio)
        << ",\n    \"maximum_one_gain_overlap_ratio\": "
        << static_cast<double>(
               report.periodic_shell_geometry
                   .maximum_one_gain_overlap_ratio)
        << ",\n    \"maximum_three_gain_overlap_ratio\": "
        << static_cast<double>(
               report.periodic_shell_geometry
                   .maximum_three_gain_overlap_ratio)
        << ",\n    \"all_bounds_hold\": "
        << (report.periodic_shell_geometry.all_bounds_hold
                ? "true"
                : "false")
        << "\n  },\n"
        << "  \"periodic_tail_bound\": {\n"
        << "    \"cutoff\": " << report.periodic_tail_bound.cutoff
        << ",\n    \"minimum_gap\": "
        << report.periodic_tail_bound.minimum_gap
        << ",\n    \"samples\": " << report.periodic_tail_bound.samples
        << ",\n    \"seed\": " << report.periodic_tail_bound.seed
        << ",\n    \"maximum_bound_ratio\": "
        << static_cast<double>(
               report.periodic_tail_bound.maximum_bound_ratio)
        << ",\n    \"nonzero_tail_seen\": "
        << (report.periodic_tail_bound.nonzero_tail_seen ? "true" : "false")
        << ",\n    \"all_bounds_hold\": "
        << (report.periodic_tail_bound.all_bounds_hold ? "true" : "false")
        << "\n  },\n"
        << "  \"far_tail_closure\": {\n"
        << "    \"base_gap\": " << report.far_tail_closure.base_gap
        << ",\n    \"samples\": " << report.far_tail_closure.samples
        << ",\n    \"seed\": " << report.far_tail_closure.seed
        << ",\n    \"viscosity\": "
        << static_cast<double>(report.far_tail_closure.viscosity)
        << ",\n    \"maximum_normalized_remainder_ratio\": "
        << static_cast<double>(
               report.far_tail_closure.maximum_normalized_remainder_ratio)
        << ",\n    \"all_bounds_hold\": "
        << (report.far_tail_closure.all_bounds_hold ? "true" : "false")
        << "\n  },\n"
        << "  \"transition_block_scaling\": {\n"
        << "    \"band_count_logarithm_power\": \""
        << report.transition_block_scaling.band_count_logarithm_power.str()
        << "\",\n    \"post_young_logarithm_power\": \""
        << report.transition_block_scaling.post_young_logarithm_power.str()
        << "\",\n    \"post_young_enstrophy_power\": \""
        << report.transition_block_scaling.post_young_enstrophy_power.str()
        << "\",\n    \"energy_time_integrable_enstrophy_power\": \""
        << report.transition_block_scaling
               .energy_time_integrable_enstrophy_power.str()
        << "\",\n    \"required_pointwise_depletion_power\": \""
        << report.transition_block_scaling
               .required_pointwise_depletion_power.str()
        << "\",\n    \"logarithmic_band_count_changes_polynomial_power\": "
        << (report.transition_block_scaling
                    .logarithmic_band_count_changes_polynomial_power
                ? "true"
                : "false")
        << ",\n    \"energy_identity_closes_transition_block\": "
        << (report.transition_block_scaling
                    .energy_identity_closes_transition_block
                ? "true"
                : "false")
        << "\n  },\n"
        << "  \"helical_triad_certificate\": {\n"
        << "    \"cutoff\": " << report.helical_triad_certificate.cutoff
        << ",\n    \"samples\": "
        << report.helical_triad_certificate.samples
        << ",\n    \"seed\": " << report.helical_triad_certificate.seed
        << ",\n    \"maximum_velocity_reconstruction_residual\": "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_velocity_reconstruction_residual)
        << ",\n    \"maximum_total_reconstruction_residual\": "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_total_reconstruction_residual)
        << ",\n    \"maximum_local_reconstruction_residual\": "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_local_reconstruction_residual)
        << ",\n    \"maximum_pure_heterochiral_absolute_local\": "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_pure_heterochiral_absolute_local)
        << ",\n    \"maximum_pure_homochiral_local_stretching\": "
        << static_cast<double>(report.helical_triad_certificate
                                   .maximum_pure_homochiral_local_stretching)
        << ",\n    \"nonzero_pure_homochiral_local_seen\": "
        << (report.helical_triad_certificate
                    .nonzero_pure_homochiral_local_seen
                ? "true"
                : "false")
        << ",\n    \"all_reconstruction_checks_hold\": "
        << (report.helical_triad_certificate
                    .all_reconstruction_checks_hold
                ? "true"
                : "false")
        << "\n  },\n"
        << "  \"fourier_galerkin\": {\n"
        << "    \"cutoff\": " << report.triad_cutoff << ",\n"
        << "    \"modes\": " << report.triad_modes << ",\n"
        << "    \"samples\": " << report.triad_samples << ",\n"
        << "    \"seed\": " << report.seed << ",\n"
        << "    \"max_normalized_energy_cancellation_residual\": "
        << static_cast<double>(report.energy_residual) << ",\n"
        << "    \"max_divergence_residual\": "
        << static_cast<double>(report.divergence_residual) << ",\n"
        << "    \"max_reality_residual\": " << static_cast<double>(report.reality_residual)
        << ",\n    \"max_classical_vortex_stretching_ratio\": "
        << static_cast<double>(report.classical_ratio)
        << ",\n    \"max_detailed_triad_cancellation_residual\": "
        << static_cast<double>(report.detailed_triad_residual)
        << ",\n    \"max_relative_detailed_triad_residual\": "
        << static_cast<double>(report.relative_detailed_triad_residual)
        << ",\n    \"locality_ratio\": 2,\n"
        << "    \"max_nonlocal_absolute_transfer_fraction\": "
        << static_cast<double>(report.nonlocal_absolute_fraction)
        << ",\n    \"max_net_flux_to_absolute_interaction_ratio\": "
        << static_cast<double>(report.flux_efficiency)
        << ",\n    \"max_local_cumulative_flux\": "
        << static_cast<double>(report.local_cumulative_flux)
        << ",\n    \"max_nonlocal_cumulative_flux\": "
        << static_cast<double>(report.nonlocal_cumulative_flux)
        << ",\n    \"max_flux_partition_residual\": "
        << static_cast<double>(report.flux_partition_residual)
        << ",\n    \"nonzero_vortex_stretching_seen\": "
        << (report.nonzero_vortex_stretching ? "true" : "false") << "\n  },\n"
        << "  \"conclusion\": \"E-Z-P monomial estimates alone cannot close the global enstrophy inequality; L4 needs an additional dynamic or geometric quantity\"\n"
        << "}\n";
}

}  // namespace lemma
