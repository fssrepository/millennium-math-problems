#include "shifted_critical_density_cli.hpp"

#include <iomanip>
#include <ostream>

namespace lemma {

void ShiftedCriticalDensityReporter::write_console(
    const ShiftedCriticalDensityReport& report, std::ostream& out) {
    const auto& value = report.diagnostic;
    out << std::setprecision(15)
        << "state=" << report.state_path
        << " cutoff=" << report.cutoff
        << " modes=" << report.modes
        << " nu=" << static_cast<double>(report.viscosity) << '\n'
        << "C_local="
        << static_cast<double>(value.local_critical_density)
        << " E0P0=" << static_cast<double>(value.initial_ep_shift)
        << " dC_dt=" << static_cast<double>(value.density_derivative)
        << " dlog_shifted_dt="
        << static_cast<double>(value.shifted_log_derivative)
        << " dlog_shifted_over_k0Z="
        << static_cast<double>(
               value.normalized_shifted_log_derivative) << '\n'
        << "dC_from_S="
        << static_cast<double>(
               report.derivative_ledger.density_from_stretching.total)
        << " dC_from_Z="
        << static_cast<double>(
               report.derivative_ledger.density_from_enstrophy.total)
        << " dC_from_P="
        << static_cast<double>(
               report.derivative_ledger.density_from_palinstrophy.total)
        << " dC_nonlinear="
        << static_cast<double>(report.derivative_ledger
               .reconstructed_density_derivative.nonlinear)
        << " dC_viscous="
        << static_cast<double>(report.derivative_ledger
               .reconstructed_density_derivative.viscous)
        << " ledger_relative_error="
        << static_cast<double>(report.derivative_ledger
               .relative_reconstruction_error) << '\n'
        << "SLD_from_S="
        << static_cast<double>(
               report.derivative_budget.normalized_from_stretching)
        << " SLD_from_Z="
        << static_cast<double>(
               report.derivative_budget.normalized_from_enstrophy)
        << " SLD_from_P="
        << static_cast<double>(
               report.derivative_budget.normalized_from_palinstrophy)
        << " SLD_nonlinear="
        << static_cast<double>(
               report.derivative_budget.normalized_nonlinear)
        << " SLD_local_nonlinear="
        << static_cast<double>(
               report.derivative_budget.normalized_local_nonlinear)
        << " SLD_nonlocal_nonlinear="
        << static_cast<double>(
               report.derivative_budget.normalized_nonlocal_nonlinear)
        << " SLD_viscous="
        << static_cast<double>(
               report.derivative_budget.normalized_viscous) << '\n'
        << "local_outer=-||grad_B_local||2="
        << static_cast<double>(
               report.quartic_identity.local_outer_state_derivative)
        << " quartic_identity_error="
        << static_cast<double>(report.quartic_identity
               .local_outer_negative_square_error) << '\n'
        << "quartic_within_shell_cancellation="
        << static_cast<double>(report.quartic_shell_ledger
               .within_shell_cancellation_fraction)
        << " quartic_between_shell_cancellation="
        << static_cast<double>(report.quartic_shell_ledger
               .between_shell_cancellation_fraction)
        << " quartic_within_mode_cancellation="
        << static_cast<double>(report.quartic_shell_ledger
               .within_target_mode_cancellation_fraction)
        << " quartic_within_eigenshell_cancellation="
        << static_cast<double>(report.quartic_shell_ledger
               .between_modes_within_eigen_shell_cancellation_fraction)
        << " quartic_within_dyadic_shell_cancellation="
        << static_cast<double>(report.quartic_shell_ledger
               .between_modes_within_dyadic_shell_cancellation_fraction)
        << " shell_envelope_max_ratio="
        << static_cast<double>(report.quartic_shell_envelope
               .maximum_bound_ratio)
        << " shell_envelope_global_ZP_ratio="
        << static_cast<double>(report.quartic_shell_envelope
               .global_zp_ratio)
        << " quartic_commutator_error="
        << static_cast<double>(report.quartic_commutator
               .identity_relative_error)
        << " quartic_projected_residual_error="
        << static_cast<double>(report.quartic_projected_residual
               .normalized_reconstruction_error)
        << " quartic_reduced_error="
        << static_cast<double>(report.quartic_reduced_ledger
               .normalized_reconstruction_error) << '\n';
}

void ShiftedCriticalDensityReporter::write_json(
    const ShiftedCriticalDensityReport& report, std::ostream& out) {
    const auto& value = report.diagnostic;
    out << std::setprecision(18)
        << "{\n  \"schema\": "
           "\"navier-stokes-shifted-critical-density-v2\",\n"
        << "  \"state\": \"" << report.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"modes\": " << report.modes << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(report.viscosity) << ",\n"
        << "  \"energy\": " << static_cast<double>(value.energy) << ",\n"
        << "  \"enstrophy\": "
        << static_cast<double>(value.enstrophy) << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(value.palinstrophy) << ",\n"
        << "  \"local_critical_density\": "
        << static_cast<double>(value.local_critical_density) << ",\n"
        << "  \"initial_E0P0_shift\": "
        << static_cast<double>(value.initial_ep_shift) << ",\n"
        << "  \"density_derivative\": "
        << static_cast<double>(value.density_derivative) << ",\n"
        << "  \"shifted_log_derivative\": "
        << static_cast<double>(value.shifted_log_derivative) << ",\n"
        << "  \"initial_frequency\": "
        << static_cast<double>(value.initial_frequency) << ",\n"
        << "  \"k0Z_normalization\": "
        << static_cast<double>(value.normalization) << ",\n"
        << "  \"normalized_shifted_log_derivative\": "
        << static_cast<double>(
               value.normalized_shifted_log_derivative) << ",\n"
        << "  \"derivative_ledger\": {\n"
        << "    \"signed_stretching\": "
        << static_cast<double>(
               report.derivative_ledger.signed_stretching) << ",\n"
        << "    \"global_signed_stretching\": "
        << static_cast<double>(report.derivative_ledger
               .global_signed_stretching) << ",\n"
        << "    \"hyperpalinstrophy\": "
        << static_cast<double>(
               report.derivative_ledger.hyperpalinstrophy) << ",\n"
        << "    \"stretching_derivative\": {\"nonlinear\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_derivative.nonlinear)
        << ", \"viscous\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_derivative.viscous)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_derivative.total) << "},\n"
        << "    \"stretching_nonlinear_roles\": {\"outer_state\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_roles.outer_state)
        << ", \"advecting_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_roles.advecting_slot)
        << ", \"advected_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_roles.advected_slot)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_roles.total) << "},\n"
        << "    \"stretching_local_nonlinear_roles\": "
           "{\"outer_state\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_local_nonlinear_roles.outer_state)
        << ", \"advecting_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_local_nonlinear_roles.advecting_slot)
        << ", \"advected_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_local_nonlinear_roles.advected_slot)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_local_nonlinear_roles.total) << "},\n"
        << "    \"stretching_nonlocal_nonlinear_roles\": "
           "{\"outer_state\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlocal_nonlinear_roles.outer_state)
        << ", \"advecting_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlocal_nonlinear_roles.advecting_slot)
        << ", \"advected_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlocal_nonlinear_roles.advected_slot)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlocal_nonlinear_roles.total) << "},\n"
        << "    \"stretching_viscous_roles\": {\"outer_state\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_viscous_roles.outer_state)
        << ", \"advecting_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_viscous_roles.advecting_slot)
        << ", \"advected_slot\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_viscous_roles.advected_slot)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_viscous_roles.total) << "},\n"
        << "    \"enstrophy_derivative\": {\"nonlinear\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_derivative.nonlinear)
        << ", \"viscous\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_derivative.viscous)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_derivative.total) << "},\n"
        << "    \"palinstrophy_derivative\": {\"nonlinear\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_derivative.nonlinear)
        << ", \"viscous\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_derivative.viscous)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_derivative.total) << "},\n"
        << "    \"density_from_stretching\": {\"nonlinear\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_stretching.nonlinear)
        << ", \"viscous\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_stretching.viscous)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_stretching.total) << "},\n"
        << "    \"density_from_enstrophy\": {\"nonlinear\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_enstrophy.nonlinear)
        << ", \"viscous\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_enstrophy.viscous)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_enstrophy.total) << "},\n"
        << "    \"density_from_palinstrophy\": {\"nonlinear\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_palinstrophy.nonlinear)
        << ", \"viscous\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_palinstrophy.viscous)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .density_from_palinstrophy.total) << "},\n"
        << "    \"reconstructed_density_derivative\": "
        << static_cast<double>(report.derivative_ledger
               .reconstructed_density_derivative.total) << ",\n"
        << "    \"oracle_density_derivative\": "
        << static_cast<double>(report.derivative_ledger
               .oracle_density_derivative) << ",\n"
        << "    \"absolute_reconstruction_error\": "
        << static_cast<double>(report.derivative_ledger
               .absolute_reconstruction_error) << ",\n"
        << "    \"relative_reconstruction_error\": "
        << static_cast<double>(report.derivative_ledger
               .relative_reconstruction_error) << ",\n"
        << "    \"stretching_nonlinear_partition\": {\"local\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_partition.local)
        << ", \"nonlocal\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_partition.nonlocal)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_partition.total) << "},\n"
        << "    \"enstrophy_nonlinear_partition\": {\"local\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_nonlinear_partition.local)
        << ", \"nonlocal\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_nonlinear_partition.nonlocal)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_nonlinear_partition.total) << "},\n"
        << "    \"palinstrophy_nonlinear_partition\": {\"local\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_nonlinear_partition.local)
        << ", \"nonlocal\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_nonlinear_partition.nonlocal)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_nonlinear_partition.total) << "},\n"
        << "    \"density_nonlinear_partition\": {\"local\": "
        << static_cast<double>(report.derivative_ledger
               .density_nonlinear_partition.local)
        << ", \"nonlocal\": "
        << static_cast<double>(report.derivative_ledger
               .density_nonlinear_partition.nonlocal)
        << ", \"total\": "
        << static_cast<double>(report.derivative_ledger
               .density_nonlinear_partition.total) << "},\n"
        << "    \"nonlinear_partition_error\": "
        << static_cast<double>(report.derivative_ledger
               .nonlinear_partition_error) << ",\n"
        << "    \"stretching_viscous_advected_cancellation\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_viscous_advected_cancellation) << ",\n"
        << "    \"enstrophy_viscous_identity_error\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_viscous_identity_error) << ",\n"
        << "    \"enstrophy_nonlinear_identity_error\": "
        << static_cast<double>(report.derivative_ledger
               .enstrophy_nonlinear_identity_error) << ",\n"
        << "    \"stretching_nonlinear_role_error\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_nonlinear_role_error) << ",\n"
        << "    \"stretching_viscous_role_error\": "
        << static_cast<double>(report.derivative_ledger
               .stretching_viscous_role_error) << ",\n"
        << "    \"palinstrophy_viscous_identity_error\": "
        << static_cast<double>(report.derivative_ledger
               .palinstrophy_viscous_identity_error) << ",\n"
        << "    \"finite\": "
        << (report.derivative_ledger.finite ? "true" : "false")
        << "\n  },\n"
        << "  \"normalized_derivative_budget\": {\n"
        << "    \"density_fraction\": "
        << static_cast<double>(
               report.derivative_budget.density_fraction) << ",\n"
        << "    \"from_stretching\": "
        << static_cast<double>(report.derivative_budget
               .normalized_from_stretching) << ",\n"
        << "    \"from_enstrophy\": "
        << static_cast<double>(report.derivative_budget
               .normalized_from_enstrophy) << ",\n"
        << "    \"from_palinstrophy\": "
        << static_cast<double>(report.derivative_budget
               .normalized_from_palinstrophy) << ",\n"
        << "    \"nonlinear\": "
        << static_cast<double>(
               report.derivative_budget.normalized_nonlinear) << ",\n"
        << "    \"local_nonlinear\": "
        << static_cast<double>(report.derivative_budget
               .normalized_local_nonlinear) << ",\n"
        << "    \"nonlocal_nonlinear\": "
        << static_cast<double>(report.derivative_budget
               .normalized_nonlocal_nonlinear) << ",\n"
        << "    \"local_quartic_budget\": {\"outer_state\": "
        << static_cast<double>(report.derivative_budget
               .normalized_local_outer_state)
        << ", \"advecting_slot\": "
        << static_cast<double>(report.derivative_budget
               .normalized_local_advecting_slot)
        << ", \"advected_slot\": "
        << static_cast<double>(report.derivative_budget
               .normalized_local_advected_slot)
        << ", \"enstrophy\": "
        << static_cast<double>(report.derivative_budget
               .normalized_local_enstrophy)
        << ", \"palinstrophy\": "
        << static_cast<double>(report.derivative_budget
               .normalized_local_palinstrophy)
        << ", \"reconstructed\": "
        << static_cast<double>(report.derivative_budget
               .reconstructed_local_nonlinear)
        << ", \"relative_error\": "
        << static_cast<double>(report.derivative_budget
               .local_nonlinear_reconstruction_error) << "},\n"
        << "    \"viscous\": "
        << static_cast<double>(
               report.derivative_budget.normalized_viscous) << ",\n"
        << "    \"reconstructed\": "
        << static_cast<double>(report.derivative_budget
               .reconstructed_normalized_rate) << ",\n"
        << "    \"relative_reconstruction_error\": "
        << static_cast<double>(report.derivative_budget
               .relative_reconstruction_error) << ",\n"
        << "    \"polynomial_SLD1P\": {\"numerator\": "
        << static_cast<double>(
               report.derivative_budget.polynomial_numerator)
        << ", \"denominator\": "
        << static_cast<double>(
               report.derivative_budget.polynomial_denominator)
        << ", \"required_coefficient\": "
        << static_cast<double>(report.derivative_budget
               .polynomial_required_coefficient)
        << ", \"relative_equivalence_error\": "
        << static_cast<double>(report.derivative_budget
               .polynomial_equivalence_error) << "},\n"
        << "    \"finite\": "
        << (report.derivative_budget.finite ? "true" : "false")
        << "\n  },\n"
        << "  \"local_quartic_identities\": {\n"
        << "    \"local_advection_H1_squared\": "
        << static_cast<double>(report.quartic_identity
               .local_advection_h1_squared) << ",\n"
        << "    \"local_outer_state_derivative\": "
        << static_cast<double>(report.quartic_identity
               .local_outer_state_derivative) << ",\n"
        << "    \"local_outer_negative_square_error\": "
        << static_cast<double>(report.quartic_identity
               .local_outer_negative_square_error) << ",\n"
        << "    \"local_enstrophy_identity_error\": "
        << static_cast<double>(report.quartic_identity
               .local_enstrophy_identity_error) << ",\n"
        << "    \"nonlocal_enstrophy_identity_error\": "
        << static_cast<double>(report.quartic_identity
               .nonlocal_enstrophy_identity_error) << ",\n"
        << "    \"finite\": "
        << (report.quartic_identity.finite ? "true" : "false")
        << "\n  },\n"
        << "  \"local_quartic_shell_ledger\": {\n"
        << "    \"normalized_total\": "
        << static_cast<double>(
               report.quartic_shell_ledger.normalized_total) << ",\n"
        << "    \"positive_dyadic_shell_total\": "
        << static_cast<double>(report.quartic_shell_ledger
               .positive_dyadic_shell_total) << ",\n"
        << "    \"negative_dyadic_shell_total\": "
        << static_cast<double>(report.quartic_shell_ledger
               .negative_dyadic_shell_total) << ",\n"
        << "    \"within_shell_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .within_shell_cancellation_fraction) << ",\n"
        << "    \"within_target_mode_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .within_target_mode_cancellation_fraction) << ",\n"
        << "    \"between_modes_within_shell_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .between_modes_within_shell_cancellation_fraction) << ",\n"
        << "    \"between_modes_within_eigen_shell_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .between_modes_within_eigen_shell_cancellation_fraction)
        << ",\n"
        << "    \"between_eigen_shells_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .between_eigen_shells_cancellation_fraction) << ",\n"
        << "    \"between_modes_within_dyadic_shell_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .between_modes_within_dyadic_shell_cancellation_fraction)
        << ",\n"
        << "    \"between_dyadic_shells_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .between_dyadic_shells_cancellation_fraction) << ",\n"
        << "    \"between_shell_cancellation_fraction\": "
        << static_cast<double>(report.quartic_shell_ledger
               .between_shell_cancellation_fraction) << ",\n"
        << "    \"raw_role_reconstruction_error\": "
        << static_cast<double>(report.quartic_shell_ledger
               .raw_role_reconstruction_error) << ",\n"
        << "    \"normalized_reconstruction_error\": "
        << static_cast<double>(report.quartic_shell_ledger
               .normalized_reconstruction_error) << ",\n"
        << "    \"shells\": [\n";
    for (std::size_t index = 0;
         index < report.quartic_shell_ledger.shells.size(); ++index) {
        const LocalQuarticShellRow& row =
            report.quartic_shell_ledger.shells[index];
        out << "      {\"shell\": " << row.shell
            << ", \"modes\": " << row.modes
            << ", \"outer_state\": "
            << static_cast<double>(row.normalized_outer_state)
            << ", \"advecting_slot\": "
            << static_cast<double>(row.normalized_advecting_slot)
            << ", \"advected_slot\": "
            << static_cast<double>(row.normalized_advected_slot)
            << ", \"enstrophy\": "
            << static_cast<double>(row.normalized_enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(row.normalized_palinstrophy)
            << ", \"total\": "
            << static_cast<double>(row.normalized_total) << "}"
            << (index + 1 == report.quartic_shell_ledger.shells.size()
                    ? "\n"
                    : ",\n");
    }
    out << "    ],\n"
        << "    \"eigen_shells\": [\n";
    bool first_eigen_shell = true;
    for (const LocalQuarticEigenShellRow& row :
         report.quartic_shell_ledger.eigen_shells) {
        if (row.modes == 0) {
            continue;
        }
        if (!first_eigen_shell) {
            out << ",\n";
        }
        first_eigen_shell = false;
        out << "      {\"wave_squared\": " << row.wave_squared
            << ", \"modes\": " << row.modes
            << ", \"outer_state\": "
            << static_cast<double>(row.normalized_outer_state)
            << ", \"advecting_slot\": "
            << static_cast<double>(row.normalized_advecting_slot)
            << ", \"advected_slot\": "
            << static_cast<double>(row.normalized_advected_slot)
            << ", \"enstrophy\": "
            << static_cast<double>(row.normalized_enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(row.normalized_palinstrophy)
            << ", \"total\": "
            << static_cast<double>(row.normalized_total) << "}";
    }
    out << "\n    ],\n"
        << "    \"dyadic_shells\": [\n";
    for (std::size_t index = 0;
         index < report.quartic_shell_ledger.dyadic_shells.size();
         ++index) {
        const LocalQuarticDyadicShellRow& row =
            report.quartic_shell_ledger.dyadic_shells[index];
        out << "      {\"shell\": " << row.shell
            << ", \"modes\": " << row.modes
            << ", \"state_energy\": "
            << static_cast<double>(row.energy)
            << ", \"state_enstrophy\": "
            << static_cast<double>(row.enstrophy)
            << ", \"state_palinstrophy\": "
            << static_cast<double>(row.palinstrophy)
            << ", \"local_advection_H1_squared\": "
            << static_cast<double>(row.local_advection_h1_squared)
            << ", \"outer_state\": "
            << static_cast<double>(row.normalized_outer_state)
            << ", \"advecting_slot\": "
            << static_cast<double>(row.normalized_advecting_slot)
            << ", \"advected_slot\": "
            << static_cast<double>(row.normalized_advected_slot)
            << ", \"enstrophy\": "
            << static_cast<double>(row.normalized_enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(row.normalized_palinstrophy)
            << ", \"total\": "
            << static_cast<double>(row.normalized_total) << "}"
            << (index + 1 ==
                        report.quartic_shell_ledger.dyadic_shells.size()
                    ? "\n"
                    : ",\n");
    }
    out << "    ],\n"
        << "    \"finite\": "
        << (report.quartic_shell_ledger.finite ? "true" : "false")
        << "\n  },\n"
        << "  \"local_quartic_shell_envelope\": {\n"
        << "    \"explicit_constant\": "
        << static_cast<double>(
               report.quartic_shell_envelope.explicit_constant) << ",\n"
        << "    \"frequency_power\": "
        << report.quartic_shell_envelope.frequency_power << ",\n"
        << "    \"required_envelope_power\": "
        << report.quartic_shell_envelope.required_envelope_power << ",\n"
        << "    \"palinstrophy_power\": "
        << report.quartic_shell_envelope.palinstrophy_power << ",\n"
        << "    \"residual_frequency_gain\": "
        << report.quartic_shell_envelope.residual_frequency_gain << ",\n"
        << "    \"maximum_bound_ratio\": "
        << static_cast<double>(report.quartic_shell_envelope
               .maximum_bound_ratio) << ",\n"
        << "    \"maximum_target_frequency_ratio\": "
        << static_cast<double>(report.quartic_shell_envelope
               .maximum_target_frequency_ratio) << ",\n"
        << "    \"maximum_advected_frequency_ratio\": "
        << static_cast<double>(report.quartic_shell_envelope
               .maximum_advected_frequency_ratio) << ",\n"
        << "    \"maximum_interaction_count_ratio\": "
        << static_cast<double>(report.quartic_shell_envelope
               .maximum_interaction_count_ratio) << ",\n"
        << "    \"neighborhood_H3_moment\": "
        << static_cast<double>(report.quartic_shell_envelope
               .neighborhood_h3_moment) << ",\n"
        << "    \"neighborhood_H4_moment\": "
        << static_cast<double>(report.quartic_shell_envelope
               .neighborhood_h4_moment) << ",\n"
        << "    \"shell_product_sum\": "
        << static_cast<double>(report.quartic_shell_envelope
               .shell_product_sum) << ",\n"
        << "    \"state_enstrophy\": "
        << static_cast<double>(report.quartic_shell_envelope
               .state_enstrophy) << ",\n"
        << "    \"state_palinstrophy\": "
        << static_cast<double>(report.quartic_shell_envelope
               .state_palinstrophy) << ",\n"
        << "    \"interpolated_H3_bound\": "
        << static_cast<double>(report.quartic_shell_envelope
               .interpolated_h3_bound) << ",\n"
        << "    \"H3_overlap_constant\": "
        << static_cast<double>(report.quartic_shell_envelope
               .h3_overlap_constant) << ",\n"
        << "    \"H4_overlap_constant\": "
        << static_cast<double>(report.quartic_shell_envelope
               .h4_overlap_constant) << ",\n"
        << "    \"actual_global_local_advection_H1_squared\": "
        << static_cast<double>(report.quartic_shell_envelope
               .actual_global_local_advection_h1_squared) << ",\n"
        << "    \"global_shell_product_bound\": "
        << static_cast<double>(report.quartic_shell_envelope
               .global_shell_product_bound) << ",\n"
        << "    \"global_ZP_bound\": "
        << static_cast<double>(report.quartic_shell_envelope
               .global_zp_bound) << ",\n"
        << "    \"global_shell_product_ratio\": "
        << static_cast<double>(report.quartic_shell_envelope
               .global_shell_product_ratio) << ",\n"
        << "    \"global_ZP_ratio\": "
        << static_cast<double>(report.quartic_shell_envelope
               .global_zp_ratio) << ",\n"
        << "    \"cutoff_independent\": "
        << (report.quartic_shell_envelope.cutoff_independent
                ? "true"
                : "false")
        << ",\n    \"all_inputs_in_neighboring_shells\": "
        << (report.quartic_shell_envelope
                    .all_inputs_in_neighboring_shells
                ? "true"
                : "false")
        << ",\n    \"all_geometry_checks_hold\": "
        << (report.quartic_shell_envelope.all_geometry_checks_hold
                ? "true"
                : "false")
        << ",\n    \"global_summation_holds\": "
        << (report.quartic_shell_envelope.global_summation_holds
                ? "true"
                : "false")
        << ",\n    \"shells\": [\n";
    for (std::size_t index = 0;
         index < report.quartic_shell_envelope.shells.size(); ++index) {
        const LocalQuarticShellEnvelopeRow& row =
            report.quartic_shell_envelope.shells[index];
        out << "      {\"shell\": " << row.shell
            << ", \"radius\": " << static_cast<double>(row.radius)
            << ", \"neighborhood_energy\": "
            << static_cast<double>(row.neighborhood_energy)
            << ", \"actual_local_advection_H1_squared\": "
            << static_cast<double>(
                   row.actual_local_advection_h1_squared)
            << ", \"explicit_bound\": "
            << static_cast<double>(row.explicit_bound)
            << ", \"bound_ratio\": "
            << static_cast<double>(row.bound_ratio) << "}"
            << (index + 1 ==
                        report.quartic_shell_envelope.shells.size()
                    ? "\n"
                    : ",\n");
    }
    out << "    ],\n"
        << "    \"all_bounds_hold\": "
        << (report.quartic_shell_envelope.all_bounds_hold
                ? "true"
                : "false")
        << "\n  },\n"
        << "  \"local_quartic_commutator\": {\n"
        << "    \"outer_state_derivative\": "
        << static_cast<double>(report.quartic_commutator
               .outer_state_derivative) << ",\n"
        << "    \"advected_slot_derivative\": "
        << static_cast<double>(report.quartic_commutator
               .advected_slot_derivative) << ",\n"
        << "    \"combined_derivative\": "
        << static_cast<double>(report.quartic_commutator
               .combined_derivative) << ",\n"
        << "    \"negative_commutator_pairing\": "
        << static_cast<double>(report.quartic_commutator
               .negative_commutator_pairing) << ",\n"
        << "    \"local_advection_L2_squared\": "
        << static_cast<double>(report.quartic_commutator
               .local_advection_l2_squared) << ",\n"
        << "    \"commutator_L2_squared\": "
        << static_cast<double>(report.quartic_commutator
               .commutator_l2_squared) << ",\n"
        << "    \"cauchy_ratio\": "
        << static_cast<double>(
               report.quartic_commutator.cauchy_ratio) << ",\n"
        << "    \"identity_relative_error\": "
        << static_cast<double>(report.quartic_commutator
               .identity_relative_error) << ",\n"
        << "    \"maximum_symbol_ratio\": "
        << static_cast<double>(report.quartic_commutator
               .maximum_symbol_ratio) << ",\n"
        << "    \"symbol_bound_holds\": "
        << (report.quartic_commutator.symbol_bound_holds
                ? "true"
                : "false")
        << ",\n    \"finite\": "
        << (report.quartic_commutator.finite ? "true" : "false")
        << "\n  },\n"
        << "  \"local_quartic_projected_residual\": {\n"
        << "    \"projected_pairing\": "
        << static_cast<double>(report.quartic_projected_residual
               .projected_pairing) << ",\n"
        << "    \"expanded_negative_square\": "
        << static_cast<double>(report.quartic_projected_residual
               .expanded_negative_square) << ",\n"
        << "    \"expanded_enstrophy_remainder\": "
        << static_cast<double>(report.quartic_projected_residual
               .expanded_enstrophy_remainder) << ",\n"
        << "    \"expanded_palinstrophy_cross\": "
        << static_cast<double>(report.quartic_projected_residual
               .expanded_palinstrophy_cross) << ",\n"
        << "    \"expanded_total\": "
        << static_cast<double>(report.quartic_projected_residual
               .expanded_total) << ",\n"
        << "    \"completed_negative_square\": "
        << static_cast<double>(report.quartic_projected_residual
               .completed_negative_square) << ",\n"
        << "    \"completed_enstrophy_remainder\": "
        << static_cast<double>(report.quartic_projected_residual
               .completed_enstrophy_remainder) << ",\n"
        << "    \"completed_hyperpalinstrophy_remainder\": "
        << static_cast<double>(report.quartic_projected_residual
               .completed_hyperpalinstrophy_remainder) << ",\n"
        << "    \"completed_total\": "
        << static_cast<double>(report.quartic_projected_residual
               .completed_total) << ",\n"
        << "    \"completion_coefficient\": "
        << static_cast<double>(report.quartic_projected_residual
               .completion_coefficient) << ",\n"
        << "    \"expansion_relative_error\": "
        << static_cast<double>(report.quartic_projected_residual
               .expansion_relative_error) << ",\n"
        << "    \"completion_relative_error\": "
        << static_cast<double>(report.quartic_projected_residual
               .completion_relative_error) << ",\n"
        << "    \"normalized_projected_pairing\": "
        << static_cast<double>(report.quartic_projected_residual
               .normalized_projected_pairing) << ",\n"
        << "    \"normalized_advecting_slot\": "
        << static_cast<double>(report.quartic_projected_residual
               .normalized_advecting_slot) << ",\n"
        << "    \"normalized_advected_slot\": "
        << static_cast<double>(report.quartic_projected_residual
               .normalized_advected_slot) << ",\n"
        << "    \"normalized_total\": "
        << static_cast<double>(report.quartic_projected_residual
               .normalized_total) << ",\n"
        << "    \"expected_normalized_local_quartet\": "
        << static_cast<double>(report.quartic_projected_residual
               .expected_normalized_local_quartet) << ",\n"
        << "    \"normalized_reconstruction_error\": "
        << static_cast<double>(report.quartic_projected_residual
               .normalized_reconstruction_error) << ",\n"
        << "    \"finite\": "
        << (report.quartic_projected_residual.finite
                ? "true"
                : "false")
        << "\n  },\n"
        << "  \"local_quartic_reduced_ledger\": {\n"
        << "    \"negative_commutator_pairing\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .negative_commutator_pairing) << ",\n"
        << "    \"enstrophy_normalization_remainder\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .enstrophy_normalization_remainder) << ",\n"
        << "    \"palinstrophy_normalization_remainder\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .palinstrophy_normalization_remainder) << ",\n"
        << "    \"reduced_pairing\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .reduced_pairing) << ",\n"
        << "    \"projected_plus_advected_pairing\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .projected_plus_advected_pairing) << ",\n"
        << "    \"raw_reconstruction_error\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .raw_reconstruction_error) << ",\n"
        << "    \"normalized_reduced_pairing\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .normalized_reduced_pairing) << ",\n"
        << "    \"normalized_advecting_slot\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .normalized_advecting_slot) << ",\n"
        << "    \"normalized_total\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .normalized_total) << ",\n"
        << "    \"expected_normalized_local_quartet\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .expected_normalized_local_quartet) << ",\n"
        << "    \"normalized_reconstruction_error\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .normalized_reconstruction_error) << ",\n"
        << "    \"polynomial_local_numerator\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .polynomial_local_numerator) << ",\n"
        << "    \"expected_polynomial_local_numerator\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .expected_polynomial_local_numerator) << ",\n"
        << "    \"polynomial_reconstruction_error\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .polynomial_reconstruction_error) << ",\n"
        << "    \"cancellation_fraction\": "
        << static_cast<double>(report.quartic_reduced_ledger
               .cancellation_fraction) << ",\n"
        << "    \"finite\": "
        << (report.quartic_reduced_ledger.finite ? "true" : "false")
        << "\n  },\n"
        << "  \"finite\": " << (value.finite ? "true" : "false")
        << "\n}\n";
}

}  // namespace lemma
