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
               value.normalized_shifted_log_derivative) << '\n';
}

void ShiftedCriticalDensityReporter::write_json(
    const ShiftedCriticalDensityReport& report, std::ostream& out) {
    const auto& value = report.diagnostic;
    out << std::setprecision(18)
        << "{\n  \"schema\": "
           "\"navier-stokes-shifted-critical-density-v1\",\n"
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
    return report.diagnostic.finite ? 0 : 2;
}

}  // namespace lemma
