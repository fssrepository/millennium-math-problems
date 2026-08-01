#include "local_signature_factor_adversary.hpp"

#include "local_signature_density.hpp"
#include "local_signature_state_factory.hpp"
#include "parallel_executor.hpp"
#include "spectral_dynamics.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

constexpr std::array<LocalSignatureStateProfile, 3> profiles{
    LocalSignatureStateProfile::decaying,
    LocalSignatureStateProfile::flat,
    LocalSignatureStateProfile::outer_half_flat};

std::uint64_t splitmix64(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

struct FactorSample {
    SpectralReal amplification_log_derivative = 0.0L;
    SpectralReal square_density_log_derivative = 0.0L;
    SpectralReal critical_log_derivative = 0.0L;
    SpectralReal derivative_product = 0.0L;
    SpectralReal simultaneous_growth_rate = 0.0L;
    std::uint64_t seed = 0;
    bool valid = false;
};

FactorSample evaluate_sample(
    int cutoff, LocalSignatureStateProfile profile, std::uint64_t seed,
    SpectralReal viscosity, SpectralReal time_step) {
    SpectralState state = LocalSignatureStateFactory::make(
        cutoff, profile, seed);
    const LocalSignatureDensitySample before =
        LocalSignatureDensity::evaluate(state);
    SpectralGalerkin galerkin;
    galerkin.configure("fft", 1);
    const SpectralDynamics dynamics(galerkin);
    dynamics.rk4_step(state, viscosity, time_step);
    const LocalSignatureDensitySample after =
        LocalSignatureDensity::evaluate(state);
    FactorSample result;
    result.seed = seed;
    if (!(before.amplification_fourth > 0.0L) ||
        !(after.amplification_fourth > 0.0L) ||
        !(before.square_signature_density > 0.0L) ||
        !(after.square_signature_density > 0.0L) ||
        !(before.critical_density > 0.0L) ||
        !(after.critical_density > 0.0L)) {
        return result;
    }
    result.amplification_log_derivative = std::log(
        after.amplification_fourth / before.amplification_fourth) /
        time_step;
    result.square_density_log_derivative = std::log(
        after.square_signature_density /
        before.square_signature_density) / time_step;
    result.critical_log_derivative = std::log(
        after.critical_density / before.critical_density) / time_step;
    result.derivative_product = result.amplification_log_derivative *
        result.square_density_log_derivative;
    result.simultaneous_growth_rate = std::min(
        result.amplification_log_derivative,
        result.square_density_log_derivative);
    result.valid = std::isfinite(result.derivative_product) &&
        std::isfinite(result.critical_log_derivative);
    return result;
}

}  // namespace

LocalSignatureFactorAdversaryReport LocalSignatureFactorAdversary::run(
    const LocalSignatureFactorAdversaryOptions& options) {
    if (options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 8 || options.samples < 1 ||
        options.workers < 1 || !(options.viscosity > 0.0L) ||
        !(options.time_step > 0.0L)) {
        throw std::invalid_argument(
            "factor adversary requires valid cutoffs, samples, workers, viscosity, and dt");
    }
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        const SpectralState warm = LocalSignatureStateFactory::make(
            cutoff, LocalSignatureStateProfile::decaying, options.seed);
        static_cast<void>(SpectralStateOps::interactions(warm));
    }

    const std::size_t cutoff_count = static_cast<std::size_t>(
        options.maximum_cutoff - options.minimum_cutoff + 1);
    const std::size_t task_count = cutoff_count * profiles.size() *
        static_cast<std::size_t>(options.samples);
    std::vector<FactorSample> samples(task_count);
    const ParallelExecutor executor(options.workers);
    executor.for_each(task_count, [&](std::size_t task) {
        const std::size_t sample =
            task % static_cast<std::size_t>(options.samples);
        const std::size_t quotient =
            task / static_cast<std::size_t>(options.samples);
        const std::size_t profile = quotient % profiles.size();
        const int cutoff = options.minimum_cutoff +
            static_cast<int>(quotient / profiles.size());
        const std::uint64_t seed = splitmix64(
            options.seed ^ splitmix64(static_cast<std::uint64_t>(cutoff)) ^
            splitmix64(static_cast<std::uint64_t>(profile + 1)) ^
            splitmix64(static_cast<std::uint64_t>(sample + 1)));
        samples[task] = evaluate_sample(
            cutoff, profiles[profile], seed, options.viscosity,
            options.time_step);
    });

    LocalSignatureFactorAdversaryReport report;
    report.workers = executor.threads();
    report.viscosity = options.viscosity;
    report.time_step = options.time_step;
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
            LocalSignatureFactorAdversaryRow row;
            row.cutoff = cutoff;
            row.profile = LocalSignatureStateFactory::name(profiles[profile]);
            row.maximum_log_derivative_product =
                -std::numeric_limits<SpectralReal>::infinity();
            row.maximum_simultaneous_growth_rate =
                -std::numeric_limits<SpectralReal>::infinity();
            row.maximum_critical_growth_rate =
                -std::numeric_limits<SpectralReal>::infinity();
            const std::size_t base =
                (static_cast<std::size_t>(cutoff - options.minimum_cutoff) *
                     profiles.size() + profile) *
                static_cast<std::size_t>(options.samples);
            for (int sample = 0; sample < options.samples; ++sample) {
                const FactorSample& value = samples[
                    base + static_cast<std::size_t>(sample)];
                if (!value.valid) {
                    continue;
                }
                ++row.valid_samples;
                row.mean_log_derivative_product +=
                    value.derivative_product;
                if (value.derivative_product < 0.0L) {
                    ++row.opposite_direction_samples;
                }
                if (value.amplification_log_derivative > 0.0L &&
                    value.square_density_log_derivative > 0.0L) {
                    ++row.simultaneous_growth_samples;
                }
                if (value.derivative_product >
                    row.maximum_log_derivative_product) {
                    row.maximum_log_derivative_product =
                        value.derivative_product;
                    row.product_seed = value.seed;
                }
                if (value.simultaneous_growth_rate >
                    row.maximum_simultaneous_growth_rate) {
                    row.maximum_simultaneous_growth_rate =
                        value.simultaneous_growth_rate;
                    row.simultaneous_growth_seed = value.seed;
                }
                if (value.critical_log_derivative >
                    row.maximum_critical_growth_rate) {
                    row.maximum_critical_growth_rate =
                        value.critical_log_derivative;
                    row.critical_growth_seed = value.seed;
                }
            }
            if (row.valid_samples > 0) {
                row.mean_log_derivative_product /=
                    static_cast<SpectralReal>(row.valid_samples);
            }
            report.universal_opposite_direction_rejected =
                report.universal_opposite_direction_rejected ||
                row.maximum_log_derivative_product > 0.0L;
            report.universal_no_simultaneous_growth_rejected =
                report.universal_no_simultaneous_growth_rejected ||
                row.simultaneous_growth_samples > 0;
            report.rows.push_back(row);
        }
    }
    return report;
}

LocalSignatureFactorAdversaryOptions
LocalSignatureFactorAdversaryCli::parse(
    int argc, char** argv, int first) {
    LocalSignatureFactorAdversaryOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--min-cutoff") {
            options.minimum_cutoff = std::stoi(next(index, name));
        } else if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--samples") {
            options.samples = std::stoi(next(index, name));
        } else if (name == "--workers") {
            options.workers = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--state-output") {
            options.state_output_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-signature-factor option: " + name);
        }
    }
    return options;
}

void LocalSignatureFactorAdversaryCli::print_help(std::ostream& out) {
    out << "Local-signature factor adversary options:\n"
        << "  --min-cutoff K       first Fourier cutoff\n"
        << "  --max-cutoff K       last Fourier cutoff\n"
        << "  --samples N          one-step states per cutoff/profile\n"
        << "  --workers N          parallel state workers\n"
        << "  --nu X               viscosity\n"
        << "  --dt X               RK4 probe step\n"
        << "  --seed N             deterministic master seed\n"
        << "  --certificate PATH   write factor-growth JSON\n"
        << "  --state-output PATH  save the strongest simultaneous-growth state\n";
}

int LocalSignatureFactorAdversaryCli::run(
    const LocalSignatureFactorAdversaryOptions& options,
    std::ostream& out) {
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "local-signature-factor requires --certificate");
    }
    const LocalSignatureFactorAdversaryReport report =
        LocalSignatureFactorAdversary::run(options);
    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write local-signature factor certificate");
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-local-signature-factor-adversary-v1\",\n"
        << "  \"workers\": " << report.workers
        << ",\n  \"samples_per_cutoff_profile\": " << options.samples
        << ",\n  \"viscosity\": "
        << static_cast<double>(report.viscosity)
        << ",\n  \"time_step\": "
        << static_cast<double>(report.time_step)
        << ",\n  \"universal_opposite_direction_rejected\": "
        << (report.universal_opposite_direction_rejected ? "true" : "false")
        << ",\n  \"universal_no_simultaneous_growth_rejected\": "
        << (report.universal_no_simultaneous_growth_rejected
                ? "true"
                : "false")
        << ",\n  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"profile\": \"" << row.profile
            << "\", \"valid_samples\": " << row.valid_samples
            << ", \"opposite_direction_samples\": "
            << row.opposite_direction_samples
            << ", \"simultaneous_growth_samples\": "
            << row.simultaneous_growth_samples
            << ", \"mean_log_derivative_product\": "
            << static_cast<double>(row.mean_log_derivative_product)
            << ", \"maximum_log_derivative_product\": "
            << static_cast<double>(row.maximum_log_derivative_product)
            << ", \"maximum_simultaneous_growth_rate\": "
            << static_cast<double>(row.maximum_simultaneous_growth_rate)
            << ", \"maximum_critical_growth_rate\": "
            << static_cast<double>(row.maximum_critical_growth_rate)
            << ", \"product_seed\": \"" << row.product_seed
            << "\", \"simultaneous_growth_seed\": \""
            << row.simultaneous_growth_seed
            << "\", \"critical_growth_seed\": \""
            << row.critical_growth_seed << "\"}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    if (!options.state_output_path.empty()) {
        const auto winner = std::max_element(
            report.rows.begin(), report.rows.end(),
            [](const LocalSignatureFactorAdversaryRow& left,
               const LocalSignatureFactorAdversaryRow& right) {
                return left.maximum_simultaneous_growth_rate <
                       right.maximum_simultaneous_growth_rate;
            });
        if (winner != report.rows.end()) {
            const LocalSignatureStateProfile profile =
                LocalSignatureStateFactory::parse(winner->profile);
            const SpectralState state = LocalSignatureStateFactory::make(
                winner->cutoff, profile,
                winner->simultaneous_growth_seed);
            SpectralStateWriter::write_tsv(
                options.state_output_path, state,
                "local-signature simultaneous-factor-growth maximizer; cutoff=" +
                    std::to_string(winner->cutoff) + " profile=" +
                    winner->profile + " seed=" +
                    std::to_string(winner->simultaneous_growth_seed));
        }
    }

    out << "local signature factor adversary workers=" << report.workers
        << " samples/profile=" << options.samples
        << "\ncutoff,profile,valid,opposite,both_grow,mean_product,"
           "max_product,max_both_growth,max_critical_growth\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ',' << row.profile << ','
            << row.valid_samples << ','
            << row.opposite_direction_samples << ','
            << row.simultaneous_growth_samples << ','
            << static_cast<double>(row.mean_log_derivative_product) << ','
            << static_cast<double>(row.maximum_log_derivative_product) << ','
            << static_cast<double>(row.maximum_simultaneous_growth_rate)
            << ','
            << static_cast<double>(row.maximum_critical_growth_rate) << '\n';
    }
    out << "universal opposite-direction claim: "
        << (report.universal_opposite_direction_rejected
                ? "REJECTED"
                : "not rejected")
        << "\nuniversal no-simultaneous-growth claim: "
        << (report.universal_no_simultaneous_growth_rejected
                ? "REJECTED"
                : "not rejected")
        << "\nCertificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
