#include "local_signature_trajectory.hpp"

#include "local_signature_density.hpp"
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

struct LogCorrelation {
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_x2 = 0.0L;
    SpectralReal sum_y2 = 0.0L;
    SpectralReal sum_xy = 0.0L;
    int count = 0;

    void add(const LocalSignatureDensitySample& sample) {
        if (!(sample.amplification_fourth > 0.0L) ||
            !(sample.square_signature_density > 0.0L)) {
            return;
        }
        const SpectralReal x = std::log(sample.amplification_fourth);
        const SpectralReal y = std::log(sample.square_signature_density);
        sum_x += x;
        sum_y += y;
        sum_x2 += x * x;
        sum_y2 += y * y;
        sum_xy += x * y;
        ++count;
    }

    [[nodiscard]] SpectralReal coefficient() const {
        if (count < 2) {
            return 0.0L;
        }
        const SpectralReal n = static_cast<SpectralReal>(count);
        const SpectralReal covariance = n * sum_xy - sum_x * sum_y;
        const SpectralReal variance_x = n * sum_x2 - sum_x * sum_x;
        const SpectralReal variance_y = n * sum_y2 - sum_y * sum_y;
        const SpectralReal denominator = std::sqrt(
            std::max(0.0L, variance_x) * std::max(0.0L, variance_y));
        return denominator > 0.0L ? covariance / denominator : 0.0L;
    }
};

void update_peaks(LocalSignatureTrajectoryRow& row,
                  const LocalSignatureDensitySample& sample, int step) {
    if (sample.amplification > row.maximum_amplification) {
        row.maximum_amplification = sample.amplification;
        row.amplification_peak_step = step;
        row.square_density_at_amplification_peak =
            sample.square_signature_density;
        row.critical_density_at_amplification_peak =
            sample.critical_density;
    }
    if (sample.square_signature_density >
        row.maximum_square_signature_density) {
        row.maximum_square_signature_density =
            sample.square_signature_density;
        row.square_density_peak_step = step;
        row.amplification_at_square_density_peak = sample.amplification;
        row.critical_density_at_square_density_peak =
            sample.critical_density;
    }
    if (sample.critical_density > row.maximum_critical_density) {
        row.maximum_critical_density = sample.critical_density;
        row.critical_peak_step = step;
        row.amplification_at_critical_peak = sample.amplification;
        row.square_density_at_critical_peak =
            sample.square_signature_density;
    }
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
    LocalSignatureDensitySample previous =
        LocalSignatureDensity::evaluate(state);
    row.initial_amplification = previous.amplification;
    row.maximum_amplification = -1.0L;
    row.maximum_critical_density = -1.0L;
    row.maximum_square_signature_density = -1.0L;
    row.maximum_factorization_residual = previous.factorization_residual;
    LogCorrelation correlation;
    correlation.add(previous);
    update_peaks(row, previous, 0);
    for (int step = 0; step < options.steps; ++step) {
        dynamics.rk4_step(state, options.viscosity, options.time_step);
        const LocalSignatureDensitySample current =
            LocalSignatureDensity::evaluate(state);
        row.critical_integral += 0.5L * options.time_step *
            (previous.critical_density + current.critical_density);
        row.square_signature_integral += 0.5L * options.time_step *
            (previous.square_signature_density +
             current.square_signature_density);
        correlation.add(current);
        update_peaks(row, current, step + 1);
        row.maximum_factorization_residual = std::max(
            row.maximum_factorization_residual,
            current.factorization_residual);
        previous = current;
    }
    row.final_amplification = previous.amplification;
    row.log_factor_correlation = correlation.coefficient();
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
            << ", \"maximum_square_signature_density\": "
            << static_cast<double>(row.maximum_square_signature_density)
            << ", \"log_factor_correlation\": "
            << static_cast<double>(row.log_factor_correlation)
            << ", \"amplification_peak_step\": "
            << row.amplification_peak_step
            << ", \"square_density_at_amplification_peak\": "
            << static_cast<double>(
                   row.square_density_at_amplification_peak)
            << ", \"critical_density_at_amplification_peak\": "
            << static_cast<double>(
                   row.critical_density_at_amplification_peak)
            << ", \"square_density_peak_step\": "
            << row.square_density_peak_step
            << ", \"amplification_at_square_density_peak\": "
            << static_cast<double>(
                   row.amplification_at_square_density_peak)
            << ", \"critical_density_at_square_density_peak\": "
            << static_cast<double>(
                   row.critical_density_at_square_density_peak)
            << ", \"critical_peak_step\": "
            << row.critical_peak_step
            << ", \"amplification_at_critical_peak\": "
            << static_cast<double>(row.amplification_at_critical_peak)
            << ", \"square_density_at_critical_peak\": "
            << static_cast<double>(row.square_density_at_critical_peak)
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
           "square_signature_integral,log_factor_correlation,"
           "A_peak_step,R2_peak_step,critical_peak_step,"
           "factorization_residual\n";
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
            << static_cast<double>(row.log_factor_correlation) << ','
            << row.amplification_peak_step << ','
            << row.square_density_peak_step << ','
            << row.critical_peak_step << ','
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
