#include "local_signature_trajectory.hpp"

#include "local_signature_objective.hpp"
#include "parallel_executor.hpp"
#include "spectral_dynamics.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

struct DensitySample {
    SpectralReal amplification = 0.0L;
    SpectralReal critical_density = 0.0L;
    SpectralReal square_signature_density = 0.0L;
    SpectralReal factorization_residual = 0.0L;
};

DensitySample sample_density(const SpectralState& state) {
    const LocalSignatureObjectiveValue signature =
        LocalSignatureObjective::evaluate(state);
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal energy = std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
        enstrophy += wave2 * energy;
        palinstrophy += wave2 * wave2 * energy;
    }
    DensitySample result;
    result.amplification = signature.signed_amplification;
    if (!(enstrophy > 0.0L) || !(palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralReal denominator = enstrophy * palinstrophy *
        palinstrophy * palinstrophy;
    const SpectralReal transfer2 = signature.signed_local_transfer *
        signature.signed_local_transfer;
    result.critical_density = transfer2 * transfer2 / denominator;
    result.square_signature_density =
        signature.squared_signature_transfer *
        signature.squared_signature_transfer / denominator;
    const SpectralReal factored =
        std::pow(signature.signed_amplification, 4.0L) *
        result.square_signature_density;
    result.factorization_residual = std::abs(
        result.critical_density - factored) /
        std::max({std::numeric_limits<SpectralReal>::min(),
                  std::abs(result.critical_density), std::abs(factored)});
    return result;
}

LocalSignatureTrajectoryRow analyze_cutoff(
    const SpectralState& input, int cutoff,
    const LocalSignatureTrajectoryOptions& options) {
    SpectralState state = input;
    const int input_cutoff = SpectralStateOps::cutoff(input);
    if (cutoff < input_cutoff) {
        state = SpectralStateFactory::project(input, cutoff);
    } else if (cutoff > input_cutoff) {
        std::mt19937_64 generator(
            UINT64_C(0x517cc1b727220a95) ^
            static_cast<std::uint64_t>(cutoff));
        state = SpectralStateFactory::lift(input, cutoff, generator);
    }
    SpectralGalerkin galerkin;
    const int cutoff_count =
        options.maximum_cutoff - options.minimum_cutoff + 1;
    const int kernel_threads = std::max(
        1, options.workers / cutoff_count);
    galerkin.configure("fft", kernel_threads);
    SpectralDynamics dynamics(galerkin);
    LocalSignatureTrajectoryRow row;
    row.cutoff = cutoff;
    DensitySample previous = sample_density(state);
    row.initial_amplification = previous.amplification;
    row.maximum_amplification = previous.amplification;
    row.maximum_critical_density = previous.critical_density;
    row.maximum_factorization_residual = previous.factorization_residual;
    for (int step = 0; step < options.steps; ++step) {
        dynamics.rk4_step(state, options.viscosity, options.time_step);
        const DensitySample current = sample_density(state);
        row.critical_integral += 0.5L * options.time_step *
            (previous.critical_density + current.critical_density);
        row.square_signature_integral += 0.5L * options.time_step *
            (previous.square_signature_density +
             current.square_signature_density);
        row.maximum_amplification = std::max(
            row.maximum_amplification, current.amplification);
        row.maximum_critical_density = std::max(
            row.maximum_critical_density, current.critical_density);
        row.maximum_factorization_residual = std::max(
            row.maximum_factorization_residual,
            current.factorization_residual);
        previous = current;
    }
    row.final_amplification = previous.amplification;
    return row;
}

}  // namespace

LocalSignatureTrajectoryReport LocalSignatureTrajectoryAnalyzer::analyze(
    const LocalSignatureTrajectoryOptions& options) {
    if (options.state_path.empty() || options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 12 || options.steps < 1 ||
        options.time_refinement < 1 || options.workers < 1 ||
        !(options.viscosity > 0.0L) ||
        !(options.time_step > 0.0L)) {
        throw std::invalid_argument(
            "invalid local signature trajectory options");
    }
    const SpectralState input = SpectralStateReader::read_tsv(
        options.state_path);
    if (options.minimum_cutoff < SpectralStateOps::cutoff(input)) {
        throw std::invalid_argument(
            "minimum trajectory cutoff cannot be below the input cutoff");
    }
    const std::size_t cutoff_count = static_cast<std::size_t>(
        options.maximum_cutoff - options.minimum_cutoff + 1);
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        SpectralState warm = input;
        if (cutoff > SpectralStateOps::cutoff(input)) {
            std::mt19937_64 generator(0);
            warm = SpectralStateFactory::lift(input, cutoff, generator);
        }
        static_cast<void>(SpectralStateOps::interactions(warm));
    }
    LocalSignatureTrajectoryReport report;
    report.state_path = options.state_path;
    report.steps = options.steps;
    report.time_refinement = options.time_refinement;
    report.viscosity = options.viscosity;
    report.time_step = options.time_step;
    report.rows.resize(cutoff_count);
    const ParallelExecutor executor(options.workers);
    report.workers = executor.threads();
    report.kernel_threads = std::max(
        1, report.workers / static_cast<int>(cutoff_count));
    executor.for_each(cutoff_count, [&](std::size_t index) {
        const int cutoff = options.minimum_cutoff +
            static_cast<int>(index);
        report.rows[index] = analyze_cutoff(input, cutoff, options);
        if (options.time_refinement > 1) {
            LocalSignatureTrajectoryOptions refined_options = options;
            refined_options.steps *= options.time_refinement;
            refined_options.time_step /= static_cast<SpectralReal>(
                options.time_refinement);
            const LocalSignatureTrajectoryRow refined = analyze_cutoff(
                input, cutoff, refined_options);
            report.rows[index].refined_critical_integral =
                refined.critical_integral;
            report.rows[index].time_refinement_relative_difference =
                std::abs(refined.critical_integral -
                         report.rows[index].critical_integral) /
                std::max({std::numeric_limits<SpectralReal>::min(),
                          std::abs(refined.critical_integral),
                          std::abs(report.rows[index].critical_integral)});
            report.rows[index].maximum_factorization_residual = std::max(
                report.rows[index].maximum_factorization_residual,
                refined.maximum_factorization_residual);
        } else {
            report.rows[index].refined_critical_integral =
                report.rows[index].critical_integral;
        }
    });
    if (report.rows.size() >= 2) {
        const SpectralReal previous =
            report.rows[report.rows.size() - 2].critical_integral;
        const SpectralReal current = report.rows.back().critical_integral;
        report.last_relative_critical_integral_difference =
            std::abs(current - previous) /
            std::max({std::numeric_limits<SpectralReal>::min(),
                      std::abs(current), std::abs(previous)});
    }
    return report;
}

LocalSignatureTrajectoryOptions LocalSignatureTrajectoryCli::parse(
    int argc, char** argv, int first) {
    LocalSignatureTrajectoryOptions options;
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
        } else if (name == "--min-cutoff") {
            options.minimum_cutoff = std::stoi(next(index, name));
        } else if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--steps") {
            options.steps = std::stoi(next(index, name));
        } else if (name == "--time-refinement") {
            options.time_refinement = std::stoi(next(index, name));
        } else if (name == "--workers") {
            options.workers = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown local-signature-trajectory option: " + name);
        }
    }
    return options;
}

void LocalSignatureTrajectoryCli::print_help(std::ostream& out) {
    out << "Local-signature trajectory options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --min-cutoff K        first Galerkin cutoff\n"
        << "  --max-cutoff K        last Galerkin cutoff\n"
        << "  --steps N             RK4 trajectory steps\n"
        << "  --time-refinement N   rerun with dt/N at fixed horizon\n"
        << "  --workers N           parallel cutoff workers\n"
        << "  --nu X                viscosity\n"
        << "  --dt X                RK4 time step\n"
        << "  --certificate PATH    write factorization JSON\n";
}

int LocalSignatureTrajectoryCli::run(
    const LocalSignatureTrajectoryOptions& options, std::ostream& out) {
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "local-signature-trajectory requires --certificate");
    }
    const LocalSignatureTrajectoryReport report =
        LocalSignatureTrajectoryAnalyzer::analyze(options);
    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write local signature trajectory certificate");
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-local-signature-trajectory-v1\",\n"
        << "  \"state\": \"" << report.state_path
        << "\",\n  \"workers\": " << report.workers
        << ",\n  \"kernel_threads_per_cutoff\": "
        << report.kernel_threads
        << ",\n  \"steps\": " << report.steps
        << ",\n  \"time_refinement\": " << report.time_refinement
        << ",\n  \"viscosity\": "
        << static_cast<double>(report.viscosity)
        << ",\n  \"time_step\": "
        << static_cast<double>(report.time_step)
        << ",\n  \"last_relative_critical_integral_difference\": "
        << static_cast<double>(
               report.last_relative_critical_integral_difference)
        << ",\n  \"exact_factorization\": "
           "\"V^4/(Z P^3) = A_sig^4 R^2/(Z P^3)\",\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"initial_amplification\": "
            << static_cast<double>(row.initial_amplification)
            << ", \"maximum_amplification\": "
            << static_cast<double>(row.maximum_amplification)
            << ", \"final_amplification\": "
            << static_cast<double>(row.final_amplification)
            << ", \"critical_integral\": "
            << static_cast<double>(row.critical_integral)
            << ", \"refined_critical_integral\": "
            << static_cast<double>(row.refined_critical_integral)
            << ", \"time_refinement_relative_difference\": "
            << static_cast<double>(
                   row.time_refinement_relative_difference)
            << ", \"square_signature_integral\": "
            << static_cast<double>(row.square_signature_integral)
            << ", \"maximum_critical_density\": "
            << static_cast<double>(row.maximum_critical_density)
            << ", \"maximum_factorization_residual\": "
            << static_cast<double>(row.maximum_factorization_residual)
            << '}'
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << "local signature trajectory workers=" << report.workers
        << " kernel_threads/cutoff=" << report.kernel_threads
        << " steps=" << report.steps
        << "\ncutoff,initial_A,max_A,final_A,critical_integral,"
           "refined_critical_integral,time_refinement_difference,"
           "square_signature_integral,factorization_residual\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ','
            << static_cast<double>(row.initial_amplification) << ','
            << static_cast<double>(row.maximum_amplification) << ','
            << static_cast<double>(row.final_amplification) << ','
            << static_cast<double>(row.critical_integral) << ','
            << static_cast<double>(row.refined_critical_integral) << ','
            << static_cast<double>(
                   row.time_refinement_relative_difference) << ','
            << static_cast<double>(row.square_signature_integral) << ','
            << static_cast<double>(row.maximum_factorization_residual)
            << '\n';
    }
    out << "last relative critical-integral difference="
        << static_cast<double>(
               report.last_relative_critical_integral_difference)
        << "\nCertificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
