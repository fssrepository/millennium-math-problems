#include "local_sld_trajectory_evaluator.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

void validate(const LocalSldTrajectoryEvaluatorOptions& options) {
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.trajectory_steps < 1 || options.threads < 1 ||
        !(options.viscosity > 0.0L) || !(options.time_step > 0.0L) ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft")) {
        throw std::invalid_argument(
            "local-sld-trajectory-evaluate requires state, certificate, and valid evolution options");
    }
}

void write_certificate(
    const LocalSldTrajectoryEvaluatorReport& report,
    const LocalSldTrajectoryEvaluatorOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write local SLD trajectory evaluation certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-trajectory-evaluation-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"backend\": \"" << options.backend << "\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(options.viscosity) << ",\n"
        << "  \"time_step\": "
        << static_cast<double>(options.time_step) << ",\n"
        << "  \"trajectory_steps\": "
        << options.trajectory_steps << ",\n"
        << "  \"initial_sld_ratio\": "
        << static_cast<double>(report.initial.terminal_ratio) << ",\n"
        << "  \"terminal_sld_ratio\": "
        << static_cast<double>(report.terminal.terminal_ratio) << ",\n"
        << "  \"maximum_sld_ratio\": "
        << static_cast<double>(report.maximum.terminal_ratio) << ",\n"
        << "  \"objective_step\": " << report.maximum.steps << ",\n"
        << "  \"objective_time\": "
        << static_cast<double>(options.time_step *
               static_cast<SpectralReal>(report.maximum.steps)) << ",\n"
        << "  \"refined_maximum_sld_ratio\": "
        << static_cast<double>(report.refined_maximum.terminal_ratio)
        << ",\n"
        << "  \"refined_objective_step\": "
        << report.refined_maximum.steps << ",\n"
        << "  \"time_step_relative_error\": "
        << static_cast<double>(report.time_step_relative_error) << ",\n"
        << "  \"initial_frequency\": "
        << static_cast<double>(report.maximum.initial_frequency) << ",\n"
        << "  \"initial_ep_shift\": "
        << static_cast<double>(report.maximum.initial_ep_shift) << ",\n"
        << "  \"finite\": " << (report.finite ? "true" : "false")
        << ",\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_evaluation_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldTrajectoryEvaluatorReport LocalSldTrajectoryEvaluator::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    SpectralReal viscosity,
    SpectralReal time_step,
    int trajectory_steps) {
    if (trajectory_steps < 1 || !(viscosity > 0.0L) ||
        !(time_step > 0.0L)) {
        throw std::invalid_argument(
            "invalid local SLD trajectory evaluation parameters");
    }
    const LocalSldTrajectoryAdjoint trajectory(dynamics);
    LocalSldTrajectoryEvaluatorReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.initial = trajectory.terminal_value(
        state, viscosity, time_step, 0);
    report.terminal = trajectory.terminal_value(
        state, viscosity, time_step, trajectory_steps);
    report.maximum = trajectory.maximum_value(
        state, viscosity, time_step, trajectory_steps);
    report.refined_maximum = trajectory.maximum_value(
        state, viscosity, 0.5L * time_step, 2 * trajectory_steps);
    report.time_step_relative_error = std::abs(
        report.refined_maximum.terminal_ratio -
        report.maximum.terminal_ratio) /
        std::max({std::abs(report.refined_maximum.terminal_ratio),
                  std::abs(report.maximum.terminal_ratio), 1e-30L});
    report.finite = report.initial.finite && report.terminal.finite &&
        report.maximum.finite && report.refined_maximum.finite &&
        std::isfinite(report.time_step_relative_error);
    return report;
}

LocalSldTrajectoryEvaluatorOptions LocalSldTrajectoryEvaluatorCli::parse(
    int argc, char** argv, int first) {
    LocalSldTrajectoryEvaluatorOptions options;
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
        } else if (name == "--backend") {
            options.backend = next(index, name);
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown local-sld-trajectory-evaluate option: " + name);
        }
    }
    validate(options);
    return options;
}

void LocalSldTrajectoryEvaluatorCli::print_help(std::ostream& out) {
    out << "Local frozen-SLD trajectory evaluator options:\n"
        << "  --state PATH         input Fourier TSV\n"
        << "  --trajectory-steps N RK4 horizon\n"
        << "  --threads N          FFT/direct workers (default 12)\n"
        << "  --nu X               viscosity (default 0.1)\n"
        << "  --dt X               RK4 step (default 0.001)\n"
        << "  --backend NAME       direct oracle, fft, or auto\n"
        << "  --certificate PATH   write English JSON evaluation\n";
}

int LocalSldTrajectoryEvaluatorCli::run(
    const LocalSldTrajectoryEvaluatorOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure(options.backend, options.threads);
    const SpectralDynamics dynamics(galerkin);
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    const LocalSldTrajectoryEvaluatorReport report =
        LocalSldTrajectoryEvaluator::evaluate(
            dynamics, state, options.viscosity, options.time_step,
            options.trajectory_steps);
    write_certificate(report, options);
    out << std::setprecision(12)
        << "local frozen-SLD trajectory cutoff=" << report.cutoff
        << " initial=" << static_cast<double>(report.initial.terminal_ratio)
        << " terminal="
        << static_cast<double>(report.terminal.terminal_ratio)
        << " maximum="
        << static_cast<double>(report.maximum.terminal_ratio)
        << " refined="
        << static_cast<double>(report.refined_maximum.terminal_ratio)
        << " peak_step=" << report.maximum.steps
        << " time_refinement="
        << static_cast<double>(report.time_step_relative_error) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.finite ? 0 : 2;
}

}  // namespace lemma
