#include "local_signature_adversary.hpp"

#include "local_triad_symmetrizer.hpp"
#include "parallel_executor.hpp"

#include <algorithm>
#include <array>
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

constexpr std::array<const char*, 3> profiles{
    "decaying", "flat", "outer-half-flat"};

std::uint64_t splitmix64(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

SpectralState make_state(int cutoff, std::size_t profile,
                         std::uint64_t seed) {
    std::mt19937_64 generator(seed);
    SpectralState state = SpectralStateFactory::random(cutoff, generator);
    if (profile == 0) {
        SpectralStateOps::normalize_energy(state);
        return state;
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        if (profile == 2 &&
            2 * std::max({std::abs(wave.x), std::abs(wave.y),
                          std::abs(wave.z)}) <= cutoff) {
            state.velocity[index] = {};
            continue;
        }
        const SpectralReal undo_decay = std::pow(
            1.0L + static_cast<SpectralReal>(norm_squared(wave)), 1.25L);
        for (SpectralComplex& component : state.velocity[index]) {
            component *= undo_decay;
        }
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

SpectralReal fitted_log_slope(
    const std::vector<LocalSignatureAdversaryRow>& rows,
    const std::string& profile, bool amplification) {
    SpectralReal count = 0.0L;
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    for (const auto& row : rows) {
        const SpectralReal value = amplification
            ? row.maximum_signed_amplification
            : row.maximum_effective_count;
        if (row.profile != profile || row.cutoff <= 1 || !(value > 0.0L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
        const SpectralReal y = std::log(value);
        count += 1.0L;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }
    const SpectralReal denominator = count * sum_xx - sum_x * sum_x;
    if (count < 2.0L || std::abs(denominator) <
            std::numeric_limits<SpectralReal>::epsilon()) {
        return 0.0L;
    }
    return (count * sum_xy - sum_x * sum_y) / denominator;
}

}  // namespace

LocalSignatureAdversaryReport LocalSignatureAdversary::run(
    const LocalSignatureAdversaryOptions& options) {
    if (options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 8 || options.samples < 1 ||
        options.workers < 1) {
        throw std::invalid_argument(
            "signature search requires 1<=min<=max<=8 and positive samples/workers");
    }

    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        const SpectralState warm = make_state(cutoff, 0, options.seed);
        static_cast<void>(SpectralStateOps::interactions(warm));
    }

    struct Sample {
        SpectralReal effective_count = 0.0L;
        SpectralReal dominant_fraction = 0.0L;
        std::size_t coherent_signatures = 0;
        SpectralReal signed_amplification = 0.0L;
        std::uint64_t seed = 0;
    };
    const std::size_t cutoff_count = static_cast<std::size_t>(
        options.maximum_cutoff - options.minimum_cutoff + 1);
    const std::size_t task_count = cutoff_count * profiles.size() *
        static_cast<std::size_t>(options.samples);
    std::vector<Sample> samples(task_count);
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
        const LocalTriadSymmetryReport report =
            LocalTriadSymmetrizer::analyze(
                make_state(cutoff, profile, seed));
        samples[task] = {report.effective_coherent_signature_count,
                         report.dominant_coherent_signature_fraction,
                         report.coherent_signature_count,
                         report.signed_signature_amplification, seed};
    });

    LocalSignatureAdversaryReport report;
    report.workers = executor.threads();
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
            LocalSignatureAdversaryRow row;
            row.cutoff = cutoff;
            row.profile = profiles[profile];
            row.samples = options.samples;
            const std::size_t base =
                (static_cast<std::size_t>(cutoff - options.minimum_cutoff) *
                     profiles.size() +
                 profile) * static_cast<std::size_t>(options.samples);
            for (int sample = 0; sample < options.samples; ++sample) {
                const Sample& value = samples[
                    base + static_cast<std::size_t>(sample)];
                row.mean_effective_count += value.effective_count;
                row.mean_signed_amplification +=
                    value.signed_amplification;
                if (value.effective_count > row.maximum_effective_count) {
                    row.maximum_effective_count = value.effective_count;
                    row.maximizer_dominant_fraction = value.dominant_fraction;
                    row.maximum_coherent_signatures =
                        value.coherent_signatures;
                    row.maximizing_seed = value.seed;
                }
                if (value.signed_amplification >
                    row.maximum_signed_amplification) {
                    row.maximum_signed_amplification =
                        value.signed_amplification;
                    row.amplifying_seed = value.seed;
                }
            }
            row.mean_effective_count /=
                static_cast<SpectralReal>(options.samples);
            row.mean_signed_amplification /=
                static_cast<SpectralReal>(options.samples);
            report.maximum_observed_effective_count = std::max(
                report.maximum_observed_effective_count,
                row.maximum_effective_count);
            report.rows.push_back(row);
        }
    }
    report.flat_maximum_log_slope = fitted_log_slope(
        report.rows, "flat", false);
    report.outer_maximum_log_slope = fitted_log_slope(
        report.rows, "outer-half-flat", false);
    report.flat_amplification_log_slope = fitted_log_slope(
        report.rows, "flat", true);
    report.outer_amplification_log_slope = fitted_log_slope(
        report.rows, "outer-half-flat", true);
    return report;
}

LocalSignatureAdversaryOptions LocalSignatureAdversaryCli::parse(
    int argc, char** argv, int first) {
    LocalSignatureAdversaryOptions options;
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
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-signature-adversary option: " + name);
        }
    }
    return options;
}

void LocalSignatureAdversaryCli::print_help(std::ostream& out) {
    out << "Local-signature adversary options:\n"
        << "  --min-cutoff K       first Fourier cutoff\n"
        << "  --max-cutoff K       last Fourier cutoff\n"
        << "  --samples N          random states per cutoff/profile\n"
        << "  --workers N          parallel state workers\n"
        << "  --seed N             deterministic master seed\n"
        << "  --certificate PATH   write search JSON\n";
}

int LocalSignatureAdversaryCli::run(
    const LocalSignatureAdversaryOptions& options, std::ostream& out) {
    if (options.certificate_path.empty()) {
        throw std::invalid_argument(
            "local-signature-adversary requires --certificate");
    }
    const LocalSignatureAdversaryReport report =
        LocalSignatureAdversary::run(options);
    const std::filesystem::path parent =
        std::filesystem::path(options.certificate_path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write local-signature adversary certificate: " +
            options.certificate_path);
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-local-signature-adversary-v1\",\n"
        << "  \"workers\": " << report.workers
        << ",\n  \"samples_per_cutoff_profile\": " << options.samples
        << ",\n  \"maximum_observed_effective_count\": "
        << static_cast<double>(report.maximum_observed_effective_count)
        << ",\n  \"flat_maximum_log_slope\": "
        << static_cast<double>(report.flat_maximum_log_slope)
        << ",\n  \"outer_maximum_log_slope\": "
        << static_cast<double>(report.outer_maximum_log_slope)
        << ",\n  \"flat_amplification_log_slope\": "
        << static_cast<double>(report.flat_amplification_log_slope)
        << ",\n  \"outer_amplification_log_slope\": "
        << static_cast<double>(report.outer_amplification_log_slope)
        << ",\n  \"finite_state_search_is_not_a_proof\": true,\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        certificate << "    {\"cutoff\": " << row.cutoff
            << ", \"profile\": \"" << row.profile
            << "\", \"samples\": " << row.samples
            << ", \"mean_effective_count\": "
            << static_cast<double>(row.mean_effective_count)
            << ", \"maximum_effective_count\": "
            << static_cast<double>(row.maximum_effective_count)
            << ", \"maximizer_dominant_fraction\": "
            << static_cast<double>(row.maximizer_dominant_fraction)
            << ", \"maximum_coherent_signatures\": "
            << row.maximum_coherent_signatures
            << ", \"maximizing_seed\": \"" << row.maximizing_seed
            << "\", \"mean_signed_amplification\": "
            << static_cast<double>(row.mean_signed_amplification)
            << ", \"maximum_signed_amplification\": "
            << static_cast<double>(row.maximum_signed_amplification)
            << ", \"amplifying_seed\": \"" << row.amplifying_seed
            << "\"}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << "local signature adversary workers=" << report.workers
        << " samples/profile=" << options.samples
        << "\ncutoff,profile,mean_N_eff,max_N_eff,dominant_fraction,"
           "coherent_signatures,seed,mean_A_sig,max_A_sig,A_seed\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ',' << row.profile << ','
            << static_cast<double>(row.mean_effective_count) << ','
            << static_cast<double>(row.maximum_effective_count) << ','
            << static_cast<double>(row.maximizer_dominant_fraction) << ','
            << row.maximum_coherent_signatures << ','
            << row.maximizing_seed << ','
            << static_cast<double>(row.mean_signed_amplification) << ','
            << static_cast<double>(row.maximum_signed_amplification) << ','
            << row.amplifying_seed << '\n';
    }
    out << "flat fitted maximum slope="
        << static_cast<double>(report.flat_maximum_log_slope)
        << " outer-half-flat slope="
        << static_cast<double>(report.outer_maximum_log_slope)
        << "\nflat amplification slope="
        << static_cast<double>(report.flat_amplification_log_slope)
        << " outer-half-flat amplification slope="
        << static_cast<double>(report.outer_amplification_log_slope) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
