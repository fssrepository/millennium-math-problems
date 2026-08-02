#include "local_sld_cyclic_trajectory_ansatz.hpp"

#include "local_sld_cyclic_basis.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numbers>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {
namespace {

void validate(const LocalSldCyclicTrajectoryOptions& options) {
    if (options.cutoff < 2 || options.cutoff > 12 ||
        options.coarse_samples < 32 ||
        options.refinement_iterations < 1 ||
        options.trajectory_steps < 1 || options.threads < 1 ||
        !(options.viscosity > 0.0L) || !(options.time_step > 0.0L) ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft")) {
        throw std::invalid_argument(
            "invalid cyclic frozen-trajectory ansatz options");
    }
}

SpectralReal projected_gradient_norm(
    SpectralIncrement gradient,
    const SpectralState& state) {
    const SpectralReal energy = SpectralStateOps::energy(state);
    const SpectralReal radial = LocalSldCyclicBasis::pairing(
        gradient, state.velocity) / energy;
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            gradient[mode][component] -=
                radial * state.velocity[mode][component];
        }
    }
    return std::sqrt(std::max(
        0.0L, LocalSldCyclicBasis::pairing(gradient, gradient)));
}

void write_certificate(
    const LocalSldCyclicTrajectoryReport& report,
    const LocalSldCyclicTrajectoryOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write cyclic trajectory ansatz certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-cyclic-trajectory-ansatz-v1\",\n"
        << "  \"ansatz\": \"normalized cyclic axis shear plus its normalized local quadratic advection response\",\n"
        << "  \"trajectory_objective\": \"maximum frozen-initial-data local SLD ratio on the sampled trajectory\",\n"
        << "  \"cutoff\": " << options.cutoff << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(options.viscosity) << ",\n"
        << "  \"time_step\": "
        << static_cast<double>(options.time_step) << ",\n"
        << "  \"trajectory_steps\": "
        << options.trajectory_steps << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"backend\": \"" << options.backend << "\",\n"
        << "  \"angle\": " << static_cast<double>(report.angle) << ",\n"
        << "  \"axis_energy_fraction\": "
        << static_cast<double>(report.axis_energy_fraction) << ",\n"
        << "  \"response_energy_fraction\": "
        << static_cast<double>(report.response_energy_fraction) << ",\n"
        << "  \"basis_inner_product\": "
        << static_cast<double>(report.basis_inner_product) << ",\n"
        << "  \"initial_static_sld_ratio\": "
        << static_cast<double>(
               report.initial_value.signed_local_sld_ratio) << ",\n"
        << "  \"maximum_sld_ratio\": "
        << static_cast<double>(report.value.terminal_ratio) << ",\n"
        << "  \"objective_step\": " << report.value.steps << ",\n"
        << "  \"objective_time\": "
        << static_cast<double>(
               options.time_step *
               static_cast<SpectralReal>(report.value.steps)) << ",\n"
        << "  \"refined_maximum_sld_ratio\": "
        << static_cast<double>(report.refined_value.terminal_ratio) << ",\n"
        << "  \"refined_objective_step\": "
        << report.refined_value.steps << ",\n"
        << "  \"time_step_relative_error\": "
        << static_cast<double>(report.time_step_relative_error) << ",\n"
        << "  \"restricted_angle_gradient\": "
        << static_cast<double>(report.restricted_gradient) << ",\n"
        << "  \"projected_full_gradient_norm\": "
        << static_cast<double>(report.projected_full_gradient_norm) << ",\n"
        << "  \"initial_frequency\": "
        << static_cast<double>(report.value.initial_frequency) << ",\n"
        << "  \"initial_ep_shift\": "
        << static_cast<double>(report.value.initial_ep_shift) << ",\n"
        << "  \"coarse_samples\": " << report.coarse_samples << ",\n"
        << "  \"refinement_iterations\": "
        << report.refinement_iterations << ",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_search_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldCyclicTrajectoryReport
LocalSldCyclicTrajectoryAnsatz::optimize(
    const LocalSldCyclicTrajectoryOptions& options) {
    validate(options);
    SpectralGalerkin galerkin;
    galerkin.configure(options.backend, 1);
    const SpectralDynamics dynamics(galerkin);
    const LocalSldTrajectoryAdjoint trajectory(dynamics);
    const LocalQuarticClosureObjective static_objective(dynamics);
    const SpectralState axis =
        LocalSldCyclicBasis::axis_state(options.cutoff);
    const SpectralState response =
        LocalSldCyclicBasis::response_state(dynamics, axis);
    const SpectralReal period =
        2.0L * std::numbers::pi_v<SpectralReal>;
    const SpectralReal spacing = period /
        static_cast<SpectralReal>(options.coarse_samples);
    std::vector<SpectralReal> values(
        static_cast<std::size_t>(options.coarse_samples));
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for schedule(static) num_threads(options.threads)
#endif
    for (int sample = 0; sample < options.coarse_samples; ++sample) {
        const SpectralReal angle = spacing *
            static_cast<SpectralReal>(sample);
        const SpectralState state = LocalSldCyclicBasis::mix(
            axis, response, angle);
        values[static_cast<std::size_t>(sample)] =
            trajectory.maximum_value(
                state, options.viscosity, options.time_step,
                options.trajectory_steps).terminal_ratio;
    }
    const auto best = std::max_element(values.begin(), values.end());
    SpectralReal best_angle = spacing * static_cast<SpectralReal>(
        std::distance(values.begin(), best));
    SpectralReal lower = best_angle - spacing;
    SpectralReal upper = best_angle + spacing;
    const SpectralReal golden =
        (std::sqrt(5.0L) - 1.0L) / 2.0L;
    auto value_at = [&](SpectralReal angle) {
        return trajectory.maximum_value(
            LocalSldCyclicBasis::mix(axis, response, angle),
            options.viscosity, options.time_step,
            options.trajectory_steps).terminal_ratio;
    };
    SpectralReal left = upper - golden * (upper - lower);
    SpectralReal right = lower + golden * (upper - lower);
    SpectralReal left_value = value_at(left);
    SpectralReal right_value = value_at(right);
    for (int iteration = 0;
         iteration < options.refinement_iterations; ++iteration) {
        if (left_value < right_value) {
            lower = left;
            left = right;
            left_value = right_value;
            right = lower + golden * (upper - lower);
            right_value = value_at(right);
        } else {
            upper = right;
            right = left;
            right_value = left_value;
            left = upper - golden * (upper - lower);
            left_value = value_at(left);
        }
    }

    LocalSldCyclicTrajectoryReport report;
    report.angle = 0.5L * (lower + upper);
    report.axis_energy_fraction =
        std::cos(report.angle) * std::cos(report.angle);
    report.response_energy_fraction =
        std::sin(report.angle) * std::sin(report.angle);
    report.basis_inner_product = LocalSldCyclicBasis::pairing(
        axis.velocity, response.velocity);
    report.state = LocalSldCyclicBasis::mix(
        axis, response, report.angle);
    report.initial_value = static_objective.evaluate(report.state);
    const QTrajectoryGradient gradient = trajectory.maximum_gradient(
        report.state, options.viscosity, options.time_step,
        options.trajectory_steps);
    report.value = trajectory.maximum_value(
        report.state, options.viscosity, options.time_step,
        options.trajectory_steps);
    report.refined_value = trajectory.maximum_value(
        report.state, options.viscosity, 0.5L * options.time_step,
        2 * options.trajectory_steps);
    const SpectralIncrement tangent =
        LocalSldCyclicBasis::angle_tangent(
            axis, response, report.angle);
    report.restricted_gradient = LocalSldCyclicBasis::pairing(
        gradient.initial_gradient, tangent);
    report.projected_full_gradient_norm = projected_gradient_norm(
        gradient.initial_gradient, report.state);
    const SpectralReal denominator = std::max(
        std::abs(report.refined_value.terminal_ratio),
        std::numeric_limits<SpectralReal>::epsilon());
    report.time_step_relative_error = std::abs(
        report.refined_value.terminal_ratio -
        report.value.terminal_ratio) / denominator;
    report.coarse_samples = options.coarse_samples;
    report.refinement_iterations = options.refinement_iterations;
    return report;
}

LocalSldCyclicTrajectoryOptions LocalSldCyclicTrajectoryCli::parse(
    int argc, char** argv, int first) {
    LocalSldCyclicTrajectoryOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--cutoff") {
            options.cutoff = std::stoi(next(index, name));
        } else if (name == "--samples") {
            options.coarse_samples = std::stoi(next(index, name));
        } else if (name == "--refinements") {
            options.refinement_iterations = std::stoi(next(index, name));
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else if (name == "--backend") {
            options.backend = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--state") {
            options.state_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-sld-trajectory-ansatz option: " + name);
        }
    }
    validate(options);
    if (options.certificate_path.empty() || options.state_path.empty()) {
        throw std::invalid_argument(
            "local-sld-trajectory-ansatz requires --certificate and --state");
    }
    return options;
}

void LocalSldCyclicTrajectoryCli::print_help(std::ostream& out) {
    out << "Local SLD cyclic frozen-trajectory ansatz options:\n"
        << "  --cutoff K           Galerkin cutoff (default 2)\n"
        << "  --samples N          coarse periodic angle samples (default 256)\n"
        << "  --refinements N      golden-section iterations (default 64)\n"
        << "  --trajectory-steps N trajectory horizon in steps (default 500)\n"
        << "  --threads N          parallel coarse-search workers (default 12)\n"
        << "  --nu VALUE           viscosity (default 0.1)\n"
        << "  --dt VALUE           RK4 time step (default 0.001)\n"
        << "  --backend NAME       direct oracle, fft, or auto (default auto)\n"
        << "  --certificate PATH   write English JSON certificate\n"
        << "  --state PATH         write optimized Fourier state\n";
}

int LocalSldCyclicTrajectoryCli::run(
    const LocalSldCyclicTrajectoryOptions& options,
    std::ostream& out) {
    const LocalSldCyclicTrajectoryReport report =
        LocalSldCyclicTrajectoryAnsatz::optimize(options);
    SpectralStateWriter::write_tsv(
        options.state_path, report.state,
        "cyclic frozen-trajectory SLD ansatz; candidate_lemma_proved=false");
    write_certificate(report, options);
    out << std::setprecision(12)
        << "cyclic trajectory ansatz cutoff=" << options.cutoff
        << " angle=" << static_cast<double>(report.angle)
        << " maximum_SLD="
        << static_cast<double>(report.value.terminal_ratio)
        << " refined_SLD="
        << static_cast<double>(report.refined_value.terminal_ratio)
        << " peak_step=" << report.value.steps
        << " restricted_gradient="
        << static_cast<double>(report.restricted_gradient)
        << " full_gradient_norm="
        << static_cast<double>(report.projected_full_gradient_norm)
        << " time_refinement="
        << static_cast<double>(report.time_step_relative_error) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
