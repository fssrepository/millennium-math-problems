#include "shifted_critical_density_cli.hpp"

#include "state_analysis.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <stdexcept>

namespace lemma {

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
    report.quartic_commutator = LocalQuarticCommutator::evaluate(
        dynamics, state, report.derivative_ledger);
    report.quartic_projected_residual =
        LocalQuarticProjectedResidual::evaluate(
            dynamics, state, report.diagnostic,
            report.derivative_ledger);
    report.quartic_reduced_ledger = LocalQuarticReducedLedger::evaluate(
        report.quartic_commutator,
        report.quartic_projected_residual,
        report.derivative_ledger);
    report.quartic_shell_ledger = LocalQuarticShellLedger::evaluate(
        dynamics, state, report.diagnostic,
        report.derivative_ledger);
    report.quartic_shell_envelope = LocalQuarticShellEnvelope::analyze(
        state, report.quartic_shell_ledger);
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
            report.quartic_commutator.finite &&
            report.quartic_projected_residual.finite &&
            report.quartic_reduced_ledger.finite &&
            report.quartic_shell_ledger.finite &&
            report.quartic_shell_envelope.all_bounds_hold &&
            report.derivative_budget.finite
        ? 0
        : 2;
}

}  // namespace lemma
