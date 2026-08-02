#include "local_sld_cyclic_krylov_ansatz.hpp"

#include "local_sld_cyclic_basis.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <array>
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

using Coefficients = std::array<SpectralReal, 3>;
using Basis = std::array<SpectralState, 3>;

void validate(const LocalSldCyclicKrylovOptions& options) {
    if (options.cutoff < 2 || options.cutoff > 12 ||
        options.warm_angle_samples < 16 || options.iterations < 1 ||
        options.line_search_steps < 1 || options.trajectory_steps < 1 ||
        options.threads < 1 || !(options.initial_step > 0.0L) ||
        !(options.viscosity > 0.0L) || !(options.time_step > 0.0L) ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft")) {
        throw std::invalid_argument("invalid cyclic Krylov ansatz options");
    }
}

SpectralReal coefficient_norm(const Coefficients& coefficients) {
    SpectralReal result = 0.0L;
    for (const SpectralReal coefficient : coefficients) {
        result += coefficient * coefficient;
    }
    return std::sqrt(result);
}

void normalize(Coefficients& coefficients) {
    const SpectralReal norm = coefficient_norm(coefficients);
    if (!(norm > 0.0L)) {
        throw std::runtime_error("zero cyclic Krylov coefficient vector");
    }
    for (SpectralReal& coefficient : coefficients) {
        coefficient /= norm;
    }
}

SpectralState state_from_coefficients(
    const Basis& basis,
    const Coefficients& coefficients) {
    SpectralState state = basis[0];
    for (ComplexVector& value : state.velocity) {
        value = {};
    }
    for (std::size_t basis_index = 0;
         basis_index < basis.size(); ++basis_index) {
        for (std::size_t mode = 0; mode < state.velocity.size(); ++mode) {
            for (std::size_t component = 0; component < 3; ++component) {
                state.velocity[mode][component] +=
                    coefficients[basis_index] *
                    basis[basis_index].velocity[mode][component];
            }
        }
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

SpectralReal maximum_gram_error(const Basis& basis) {
    SpectralReal error = 0.0L;
    for (std::size_t left = 0; left < basis.size(); ++left) {
        for (std::size_t right = 0; right < basis.size(); ++right) {
            const SpectralReal expected = left == right ? 1.0L : 0.0L;
            error = std::max(error, std::abs(
                LocalSldCyclicBasis::pairing(
                    basis[left].velocity, basis[right].velocity) -
                expected));
        }
    }
    return error;
}

SpectralReal full_projected_gradient_norm(
    SpectralIncrement gradient,
    const SpectralState& state) {
    const SpectralReal radial = LocalSldCyclicBasis::pairing(
        gradient, state.velocity) / SpectralStateOps::energy(state);
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            gradient[mode][component] -=
                radial * state.velocity[mode][component];
        }
    }
    return std::sqrt(std::max(
        0.0L, LocalSldCyclicBasis::pairing(gradient, gradient)));
}

Coefficients restricted_gradient(
    const Basis& basis,
    const Coefficients& coefficients,
    const SpectralIncrement& state_gradient) {
    Coefficients result{};
    SpectralReal radial = 0.0L;
    for (std::size_t index = 0; index < basis.size(); ++index) {
        result[index] = LocalSldCyclicBasis::pairing(
            state_gradient, basis[index].velocity);
        radial += result[index] * coefficients[index];
    }
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] -= radial * coefficients[index];
    }
    return result;
}

void write_certificate(
    const LocalSldCyclicKrylovReport& report,
    const LocalSldCyclicKrylovOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write cyclic Krylov ansatz certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-cyclic-krylov-ansatz-v1\",\n"
        << "  \"basis\": [\"cyclic axis shear\", \"quadratic local response\", \"orthogonalized advection JVP response\"],\n"
        << "  \"cutoff\": " << options.cutoff << ",\n"
        << "  \"backend\": \"" << options.backend << "\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(options.viscosity) << ",\n"
        << "  \"time_step\": "
        << static_cast<double>(options.time_step) << ",\n"
        << "  \"trajectory_steps\": "
        << options.trajectory_steps << ",\n"
        << "  \"coefficients\": ["
        << static_cast<double>(report.coefficients[0]) << ", "
        << static_cast<double>(report.coefficients[1]) << ", "
        << static_cast<double>(report.coefficients[2]) << "],\n"
        << "  \"energy_fractions\": ["
        << static_cast<double>(report.energy_fractions[0]) << ", "
        << static_cast<double>(report.energy_fractions[1]) << ", "
        << static_cast<double>(report.energy_fractions[2]) << "],\n"
        << "  \"maximum_gram_error\": "
        << static_cast<double>(report.maximum_gram_error) << ",\n"
        << "  \"maximum_sld_ratio\": "
        << static_cast<double>(report.value.terminal_ratio) << ",\n"
        << "  \"objective_step\": " << report.value.steps << ",\n"
        << "  \"objective_time\": "
        << static_cast<double>(options.time_step *
               static_cast<SpectralReal>(report.value.steps)) << ",\n"
        << "  \"refined_maximum_sld_ratio\": "
        << static_cast<double>(report.refined_value.terminal_ratio) << ",\n"
        << "  \"refined_objective_step\": "
        << report.refined_value.steps << ",\n"
        << "  \"time_step_relative_error\": "
        << static_cast<double>(report.time_step_relative_error) << ",\n"
        << "  \"restricted_gradient_norm\": "
        << static_cast<double>(report.restricted_gradient_norm) << ",\n"
        << "  \"projected_full_gradient_norm\": "
        << static_cast<double>(report.projected_full_gradient_norm) << ",\n"
        << "  \"accepted_steps\": " << report.accepted_steps << ",\n"
        << "  \"evaluations\": " << report.evaluations << ",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_search_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldCyclicKrylovReport LocalSldCyclicKrylovAnsatz::optimize(
    const LocalSldCyclicKrylovOptions& options) {
    validate(options);
    SpectralGalerkin sampler_galerkin;
    sampler_galerkin.configure(options.backend, 1);
    const SpectralDynamics sampler_dynamics(sampler_galerkin);
    Basis basis;
    basis[0] = LocalSldCyclicBasis::axis_state(options.cutoff);
    basis[1] = LocalSldCyclicBasis::response_state(
        sampler_dynamics, basis[0]);
    basis[2] = LocalSldCyclicBasis::cubic_response_state(
        sampler_dynamics, basis[0], basis[1]);
    const LocalSldTrajectoryAdjoint sampler_trajectory(
        sampler_dynamics);
    const SpectralReal spacing =
        2.0L * std::numbers::pi_v<SpectralReal> /
        static_cast<SpectralReal>(options.warm_angle_samples);
    std::vector<SpectralReal> warm_values(
        static_cast<std::size_t>(options.warm_angle_samples));
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for schedule(static) num_threads(options.threads)
#endif
    for (int sample = 0; sample < options.warm_angle_samples; ++sample) {
        const SpectralReal angle = spacing *
            static_cast<SpectralReal>(sample);
        Coefficients trial{
            std::cos(angle), std::sin(angle), 0.0L};
        warm_values[static_cast<std::size_t>(sample)] =
            sampler_trajectory.maximum_value(
                state_from_coefficients(basis, trial),
                options.viscosity, options.time_step,
                options.trajectory_steps).terminal_ratio;
    }
    const std::size_t best_sample = static_cast<std::size_t>(
        std::distance(warm_values.begin(),
                      std::max_element(
                          warm_values.begin(), warm_values.end())));
    const SpectralReal warm_angle = spacing *
        static_cast<SpectralReal>(best_sample);
    Coefficients coefficients{
        std::cos(warm_angle), std::sin(warm_angle), 0.0L};

    SpectralGalerkin optimizer_galerkin;
    optimizer_galerkin.configure(options.backend, options.threads);
    const SpectralDynamics optimizer_dynamics(optimizer_galerkin);
    const LocalSldTrajectoryAdjoint trajectory(optimizer_dynamics);
    SpectralState state = state_from_coefficients(basis, coefficients);
    SpectralReal objective = trajectory.maximum_value(
        state, options.viscosity, options.time_step,
        options.trajectory_steps).terminal_ratio;
    SpectralReal step_size = options.initial_step;
    int accepted_steps = 0;
    int evaluations = options.warm_angle_samples + 1;
    for (int iteration = 0; iteration < options.iterations; ++iteration) {
        const QTrajectoryGradient state_gradient =
            trajectory.maximum_gradient(
                state, options.viscosity, options.time_step,
                options.trajectory_steps);
        const Coefficients gradient = restricted_gradient(
            basis, coefficients, state_gradient.initial_gradient);
        if (coefficient_norm(gradient) < 1e-13L) {
            break;
        }
        bool accepted = false;
        SpectralReal trial_step = step_size;
        for (int line = 0; line < options.line_search_steps; ++line) {
            Coefficients trial = coefficients;
            for (std::size_t index = 0; index < trial.size(); ++index) {
                trial[index] += trial_step * gradient[index];
            }
            normalize(trial);
            const SpectralState trial_state =
                state_from_coefficients(basis, trial);
            const SpectralReal trial_objective =
                trajectory.maximum_value(
                    trial_state, options.viscosity, options.time_step,
                    options.trajectory_steps).terminal_ratio;
            ++evaluations;
            if (trial_objective > objective) {
                coefficients = trial;
                state = trial_state;
                objective = trial_objective;
                step_size = std::min(1.0L, 1.5L * trial_step);
                ++accepted_steps;
                accepted = true;
                break;
            }
            trial_step *= 0.5L;
        }
        if (!accepted) {
            break;
        }
    }

    const QTrajectoryGradient final_gradient = trajectory.maximum_gradient(
        state, options.viscosity, options.time_step,
        options.trajectory_steps);
    const Coefficients coefficient_gradient = restricted_gradient(
        basis, coefficients, final_gradient.initial_gradient);
    LocalSldCyclicKrylovReport report;
    report.state = state;
    report.value = trajectory.maximum_value(
        state, options.viscosity, options.time_step,
        options.trajectory_steps);
    report.refined_value = trajectory.maximum_value(
        state, options.viscosity, 0.5L * options.time_step,
        2 * options.trajectory_steps);
    report.coefficients = coefficients;
    for (std::size_t index = 0; index < coefficients.size(); ++index) {
        report.energy_fractions[index] =
            coefficients[index] * coefficients[index];
    }
    report.maximum_gram_error = maximum_gram_error(basis);
    report.restricted_gradient_norm = coefficient_norm(
        coefficient_gradient);
    report.projected_full_gradient_norm = full_projected_gradient_norm(
        final_gradient.initial_gradient, state);
    const SpectralReal denominator = std::max(
        std::abs(report.refined_value.terminal_ratio),
        std::numeric_limits<SpectralReal>::epsilon());
    report.time_step_relative_error = std::abs(
        report.refined_value.terminal_ratio -
        report.value.terminal_ratio) / denominator;
    report.accepted_steps = accepted_steps;
    report.evaluations = evaluations;
    return report;
}

LocalSldCyclicKrylovOptions LocalSldCyclicKrylovCli::parse(
    int argc, char** argv, int first) {
    LocalSldCyclicKrylovOptions options;
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
        } else if (name == "--warm-samples") {
            options.warm_angle_samples = std::stoi(next(index, name));
        } else if (name == "--iterations") {
            options.iterations = std::stoi(next(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next(index, name));
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next(index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next(index, name));
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
                "unknown local-sld-krylov-ansatz option: " + name);
        }
    }
    validate(options);
    if (options.certificate_path.empty() || options.state_path.empty()) {
        throw std::invalid_argument(
            "local-sld-krylov-ansatz requires --certificate and --state");
    }
    return options;
}

void LocalSldCyclicKrylovCli::print_help(std::ostream& out) {
    out << "Local SLD cyclic three-vector Krylov ansatz options:\n"
        << "  --cutoff K           Galerkin cutoff (default 3)\n"
        << "  --warm-samples N     parallel two-vector angle samples\n"
        << "  --iterations N       adjoint coefficient steps\n"
        << "  --line-search N      backtracking trials per step\n"
        << "  --trajectory-steps N trajectory horizon in RK4 steps\n"
        << "  --threads N          parallel/FFT workers (default 12)\n"
        << "  --step X             initial coefficient-space step\n"
        << "  --nu X               viscosity (default 0.1)\n"
        << "  --dt X               RK4 step (default 0.001)\n"
        << "  --backend NAME       direct oracle, fft, or auto\n"
        << "  --certificate PATH   write English JSON certificate\n"
        << "  --state PATH         write optimized Fourier state\n";
}

int LocalSldCyclicKrylovCli::run(
    const LocalSldCyclicKrylovOptions& options,
    std::ostream& out) {
    const LocalSldCyclicKrylovReport report =
        LocalSldCyclicKrylovAnsatz::optimize(options);
    SpectralStateWriter::write_tsv(
        options.state_path, report.state,
        "cyclic three-vector Krylov trajectory ansatz; candidate_lemma_proved=false");
    write_certificate(report, options);
    out << std::setprecision(12)
        << "cyclic Krylov ansatz cutoff=" << options.cutoff
        << " maximum_SLD="
        << static_cast<double>(report.value.terminal_ratio)
        << " refined_SLD="
        << static_cast<double>(report.refined_value.terminal_ratio)
        << " peak_step=" << report.value.steps
        << " coefficients="
        << static_cast<double>(report.coefficients[0]) << ','
        << static_cast<double>(report.coefficients[1]) << ','
        << static_cast<double>(report.coefficients[2])
        << " restricted_gradient="
        << static_cast<double>(report.restricted_gradient_norm)
        << " full_gradient="
        << static_cast<double>(report.projected_full_gradient_norm)
        << " accepted=" << report.accepted_steps << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
