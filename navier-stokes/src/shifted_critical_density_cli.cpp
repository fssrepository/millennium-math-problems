#include "shifted_critical_density_cli.hpp"

#include "state_analysis.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

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
               .local_outer_negative_square_error) << '\n';
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
        << "  \"finite\": " << (value.finite ? "true" : "false")
        << "\n}\n";
}

ShiftedCriticalDensityOptions ShiftedCriticalDensityCli::parse(
    int argc, char** argv, int first) {
    ShiftedCriticalDensityOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--backend") {
            options.backend = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown shifted-density option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty()) {
        throw std::invalid_argument(
            "shifted-density requires --state and --certificate");
    }
    if (!(options.viscosity > 0.0L) ||
        !std::isfinite(options.viscosity) ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "shifted-density numerical options are outside their range");
    }
    if (options.backend != "auto" && options.backend != "direct" &&
        options.backend != "fft") {
        throw std::invalid_argument(
            "shifted-density backend must be auto, direct, or fft");
    }
    return options;
}

void ShiftedCriticalDensityCli::print_help(std::ostream& out) {
    out << "Shifted critical-density diagnostic options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write instantaneous derivative JSON\n"
        << "  --nu VALUE            viscosity (default 0.1)\n"
        << "  --threads N           FFT worker threads (default 12)\n"
        << "  --backend NAME        auto, direct, or fft\n";
}

int ShiftedCriticalDensityCli::run(
    const ShiftedCriticalDensityOptions& options, std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure(options.backend, options.threads);
    const SpectralDynamics dynamics(galerkin);
    const SpectralObjective objective(dynamics);
    const SpectralState state =
        SpectralStateReader::read_tsv(options.state_path);
    ShiftedCriticalDensityReport report;
    report.state_path = options.state_path;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.modes = static_cast<int>(state.waves.size());
    report.viscosity = options.viscosity;
    report.diagnostic = ShiftedCriticalDensityAnalyzer::evaluate(
        dynamics, objective, state, options.viscosity);
    report.derivative_ledger = LocalCriticalDerivativeLedger::evaluate(
        dynamics, objective, state, options.viscosity,
        TriadPartition::local, options.threads);
    report.quartic_identity = LocalQuarticIdentityLedger::evaluate(
        dynamics, state, report.derivative_ledger);
    report.derivative_budget = ShiftedCriticalDensityBudgetAnalyzer::evaluate(
        report.diagnostic, report.derivative_ledger);
    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write shifted-density certificate");
    }
    ShiftedCriticalDensityReporter::write_json(report, certificate);
    ShiftedCriticalDensityReporter::write_console(report, out);
    out << "Certificate written to " << options.certificate_path << '\n';
    return report.diagnostic.finite && report.derivative_ledger.finite &&
            report.quartic_identity.finite &&
            report.derivative_budget.finite
        ? 0
        : 2;
}

}  // namespace lemma
