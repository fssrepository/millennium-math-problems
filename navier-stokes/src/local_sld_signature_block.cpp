#include "local_sld_signature_block.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

SpectralReal pairing(const SpectralIncrement& left,
                     const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "signature block increment layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralReal relative_error(SpectralReal difference,
                            SpectralReal scale) {
    return std::abs(difference) /
        std::max(std::abs(scale), 1e-30L);
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source) {
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal weight = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

SpectralReal reconstruction_error(
    const SpectralIncrement& full,
    const SpectralIncrement& dominant,
    const SpectralIncrement& remainder) {
    SpectralReal residual2 = 0.0L;
    SpectralReal scale2 = 0.0L;
    for (std::size_t mode = 0; mode < full.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            residual2 += std::norm(
                full[mode][component] - dominant[mode][component] -
                remainder[mode][component]);
            scale2 += std::norm(full[mode][component]);
        }
    }
    return std::sqrt(residual2 /
        std::max(scale2, std::numeric_limits<SpectralReal>::min()));
}

std::array<SpectralInteger, 3> parse_signature(
    const std::string& text) {
    std::array<SpectralInteger, 3> result{};
    std::stringstream input(text);
    std::string token;
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (!std::getline(input, token, ',')) {
            throw std::invalid_argument(
                "signature must contain three comma-separated integers");
        }
        result[index] = std::stoll(token);
        if (result[index] <= 0) {
            throw std::invalid_argument(
                "signature squared lengths must be positive");
        }
    }
    if (std::getline(input, token, ',')) {
        throw std::invalid_argument(
            "signature must contain exactly three integers");
    }
    std::sort(result.begin(), result.end());
    return result;
}

void write_certificate(
    const LocalSldSignatureBlockReport& report,
    const LocalSldSignatureBlockCliOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot write local SLD signature-block certificate");
    }
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-signature-block-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"backend\": \"" << options.backend << "\",\n"
        << "  \"squared_length_signature\": ["
        << report.squared_lengths[0] << ", "
        << report.squared_lengths[1] << ", "
        << report.squared_lengths[2] << "],\n"
        << "  \"equal_low_doubling_family\": "
        << (report.equal_low_doubling_family ? "true" : "false")
        << ",\n"
        << "  \"frozen_initial_normalization\": "
        << (report.frozen_initial_normalization ? "true" : "false")
        << ",\n"
        << "  \"normalization_initial_frequency\": "
        << static_cast<double>(report.normalization_initial_frequency)
        << ",\n"
        << "  \"normalization_initial_ep_shift\": "
        << static_cast<double>(report.normalization_initial_ep_shift)
        << ",\n"
        << "  \"evolved_steps\": " << report.evolved_steps << ",\n"
        << "  \"evolved_time\": "
        << static_cast<double>(report.evolved_time) << ",\n"
        << "  \"viscosity\": "
        << static_cast<double>(report.viscosity) << ",\n"
        << "  \"dominant_interactions\": "
        << report.dominant_interactions << ",\n"
        << "  \"remainder_interactions\": "
        << report.remainder_interactions << ",\n"
        << "  \"advection_reconstruction_error\": "
        << static_cast<double>(report.advection_reconstruction_error)
        << ",\n"
        << "  \"stretching_reconstruction_error\": "
        << static_cast<double>(report.stretching_reconstruction_error)
        << ",\n"
        << "  \"bracket_reconstruction_error\": "
        << static_cast<double>(report.bracket_reconstruction_error)
        << ",\n"
        << "  \"quotient_reconstruction_error\": "
        << static_cast<double>(report.quotient_reconstruction_error)
        << ",\n"
        << "  \"dominant_advection_norm2\": "
        << static_cast<double>(report.dominant_advection_norm2) << ",\n"
        << "  \"remainder_advection_norm2\": "
        << static_cast<double>(report.remainder_advection_norm2) << ",\n"
        << "  \"advection_cross_pairing\": "
        << static_cast<double>(report.advection_cross_pairing) << ",\n"
        << "  \"full_stretching\": "
        << static_cast<double>(report.full_stretching) << ",\n"
        << "  \"dominant_stretching\": "
        << static_cast<double>(report.dominant_stretching) << ",\n"
        << "  \"remainder_stretching\": "
        << static_cast<double>(report.remainder_stretching) << ",\n"
        << "  \"full_bracket\": "
        << static_cast<double>(report.full_bracket) << ",\n"
        << "  \"dominant_closed_bracket\": "
        << static_cast<double>(report.dominant_closed_bracket) << ",\n"
        << "  \"remainder_closed_bracket\": "
        << static_cast<double>(report.remainder_closed_bracket) << ",\n"
        << "  \"cross_bracket\": "
        << static_cast<double>(report.cross_bracket) << ",\n"
        << "  \"full_signed_sld_ratio\": "
        << static_cast<double>(report.full_signed_sld_ratio) << ",\n"
        << "  \"dominant_closed_sld_ratio\": "
        << static_cast<double>(report.dominant_closed_sld_ratio) << ",\n"
        << "  \"remainder_closed_sld_ratio\": "
        << static_cast<double>(report.remainder_closed_sld_ratio) << ",\n"
        << "  \"cross_sld_ratio\": "
        << static_cast<double>(report.cross_sld_ratio) << ",\n"
        << "  \"dominant_absolute_fraction\": "
        << static_cast<double>(report.dominant_absolute_fraction) << ",\n"
        << "  \"exact_decomposition\": "
        << (report.exact_decomposition ? "true" : "false") << ",\n"
        << "  \"candidate_block_bound_proved\": false,\n"
        << "  \"finite_measurement_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldSignatureBlockReport LocalSldSignatureBlock::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    std::array<SpectralInteger, 3> squared_lengths,
    bool equal_low_doubling_family,
    SpectralReal initial_frequency,
    SpectralReal initial_ep_shift) {
    std::sort(squared_lengths.begin(), squared_lengths.end());
    const TriadSelection dominant_selection = equal_low_doubling_family
        ? TriadSelection::local_equal_low_doubling()
        : TriadSelection::local_signature(
              squared_lengths[0], squared_lengths[1], squared_lengths[2]);
    const TriadSelection remainder_selection = equal_low_doubling_family
        ? TriadSelection::local_without_equal_low_doubling()
        : TriadSelection::local_without_signature(
              squared_lengths[0], squared_lengths[1], squared_lengths[2]);
    const LocalQuarticClosureObjective full_objective(
        dynamics, TriadPartition::local);
    const LocalQuarticClosureObjective dominant_objective(
        dynamics, dominant_selection);
    const LocalQuarticClosureObjective remainder_objective(
        dynamics, remainder_selection);

    LocalSldSignatureBlockReport report;
    report.squared_lengths = squared_lengths;
    report.equal_low_doubling_family = equal_low_doubling_family;
    report.full = full_objective.evaluate(state);
    report.dominant = dominant_objective.evaluate(state);
    report.remainder = remainder_objective.evaluate(state);
    const SpectralIncrement full_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::local);
    const SpectralIncrement dominant_advection =
        dynamics.advection_direct_partition(state, dominant_selection);
    const SpectralIncrement remainder_advection =
        dynamics.advection_direct_partition(state, remainder_selection);
    report.advection_reconstruction_error = reconstruction_error(
        full_advection, dominant_advection, remainder_advection);
    report.dominant_advection_norm2 = pairing(
        dominant_advection, dominant_advection);
    report.remainder_advection_norm2 = pairing(
        remainder_advection, remainder_advection);
    report.advection_cross_pairing = pairing(
        dominant_advection, remainder_advection);

    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (TriadPartitioner::includes(
                state, interaction, dominant_selection)) {
            ++report.dominant_interactions;
        } else if (TriadPartitioner::includes(
                       state, interaction, remainder_selection)) {
            ++report.remainder_interactions;
        }
    }
    report.full_stretching = report.full.signed_stretching;
    report.dominant_stretching = report.dominant.signed_stretching;
    report.remainder_stretching = report.remainder.signed_stretching;
    report.stretching_reconstruction_error = relative_error(
        report.full_stretching - report.dominant_stretching -
            report.remainder_stretching,
        report.full_stretching);
    report.full_bracket = report.full.signed_two_entry_bracket;
    report.dominant_closed_bracket =
        report.dominant.signed_two_entry_bracket;
    report.remainder_closed_bracket =
        report.remainder.signed_two_entry_bracket;
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity);
    const SpectralIncrement dominant_ab = laplacian_weight(
        state, dominant_advection);
    const SpectralIncrement remainder_ab = laplacian_weight(
        state, remainder_advection);
    const SpectralIncrement dominant_transported_au =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, au, dominant_selection);
    const SpectralIncrement remainder_transported_au =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, au, remainder_selection);
    const SpectralIncrement dominant_remainder_advecting =
        dynamics.advection_bilinear_direct_partition(
            state, remainder_advection, state.velocity,
            dominant_selection);
    const SpectralIncrement remainder_dominant_advecting =
        dynamics.advection_bilinear_direct_partition(
            state, dominant_advection, state.velocity,
            remainder_selection);
    const SpectralReal dominant_cross = pairing(
        dominant_ab, au);
    const SpectralReal remainder_cross = pairing(
        remainder_ab, au);
    report.cross_bracket =
        -pairing(dominant_advection, remainder_ab) -
        pairing(remainder_advection, dominant_ab) +
        pairing(dominant_advection, remainder_transported_au) +
        pairing(remainder_advection, dominant_transported_au) +
        report.dominant_stretching * report.remainder_stretching /
            report.full.enstrophy +
        3.0L *
            (report.dominant_stretching * remainder_cross +
             report.remainder_stretching * dominant_cross) /
            (2.0L * report.full.palinstrophy) -
        pairing(au, dominant_remainder_advecting) -
        pairing(au, remainder_dominant_advecting);
    report.bracket_reconstruction_error = relative_error(
        report.full_bracket - report.dominant_closed_bracket -
            report.remainder_closed_bracket - report.cross_bracket,
        report.full_bracket);

    report.frozen_initial_normalization =
        initial_frequency > 0.0L && initial_ep_shift > 0.0L;
    report.normalization_initial_frequency =
        report.frozen_initial_normalization
        ? initial_frequency
        : report.full.initial_frequency;
    report.normalization_initial_ep_shift =
        report.frozen_initial_normalization
        ? initial_ep_shift
        : report.full.initial_ep_shift;
    const SpectralReal s2 = report.full.signed_stretching *
        report.full.signed_stretching;
    const SpectralReal s4 = s2 * s2;
    const SpectralReal p2 = report.full.palinstrophy *
        report.full.palinstrophy;
    const SpectralReal p4 = p2 * p2;
    const SpectralReal common_denominator =
        report.normalization_initial_frequency *
        (s4 * report.full.enstrophy * report.full.enstrophy *
             report.full.palinstrophy +
         report.normalization_initial_ep_shift *
             report.full.enstrophy * report.full.enstrophy *
             report.full.enstrophy * p4);
    report.full_signed_sld_ratio = common_denominator > 0.0L
        ? 4.0L * report.full.signed_stretching * s2 *
              report.full.enstrophy * report.full.palinstrophy *
              report.full_bracket / common_denominator
        : 0.0L;
    const SpectralReal common_multiplier = common_denominator > 0.0L
        ? 4.0L * report.full.signed_stretching *
              report.full.signed_stretching *
              report.full.signed_stretching *
              report.full.enstrophy * report.full.palinstrophy /
              common_denominator
        : 0.0L;
    report.dominant_closed_sld_ratio = common_multiplier *
        report.dominant_closed_bracket;
    report.remainder_closed_sld_ratio = common_multiplier *
        report.remainder_closed_bracket;
    report.cross_sld_ratio = common_multiplier * report.cross_bracket;
    const SpectralReal reconstructed_ratio =
        report.dominant_closed_sld_ratio +
        report.remainder_closed_sld_ratio + report.cross_sld_ratio;
    report.quotient_reconstruction_error = relative_error(
        report.full_signed_sld_ratio - reconstructed_ratio,
        report.full_signed_sld_ratio);
    const SpectralReal absolute_ratio_sum =
        std::abs(report.dominant_closed_sld_ratio) +
        std::abs(report.remainder_closed_sld_ratio) +
        std::abs(report.cross_sld_ratio);
    if (absolute_ratio_sum > 0.0L) {
        report.dominant_absolute_fraction =
            std::abs(report.dominant_closed_sld_ratio) /
            absolute_ratio_sum;
    }
    report.exact_decomposition =
        report.advection_reconstruction_error < 1e-14L &&
        report.stretching_reconstruction_error < 1e-14L &&
        report.bracket_reconstruction_error < 1e-14L &&
        report.quotient_reconstruction_error < 1e-14L;
    report.finite = report.full.finite && report.dominant.finite &&
        report.remainder.finite && report.exact_decomposition;
    return report;
}

LocalSldSignatureBlockCliOptions LocalSldSignatureBlockCli::parse(
    int argc, char** argv, int first) {
    LocalSldSignatureBlockCliOptions options;
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
        } else if (name == "--signature") {
            options.squared_lengths = parse_signature(next(index, name));
        } else if (name == "--doubling-family") {
            options.equal_low_doubling_family = true;
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--evolve-steps") {
            options.evolve_steps = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else if (name == "--backend") {
            options.backend = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-sld-block option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.threads < 1 || options.threads > 256 ||
        options.evolve_steps < 0 || !(options.viscosity > 0.0L) ||
        !(options.time_step > 0.0L) ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft")) {
        throw std::invalid_argument(
            "local-sld-block requires --state, --certificate, and valid threads");
    }
    return options;
}

void LocalSldSignatureBlockCli::print_help(std::ostream& out) {
    out << "Local SLD exact signature-block options:\n"
        << "  --state PATH         input Fourier TSV\n"
        << "  --signature A,B,C    squared-length signature (default 1,1,2)\n"
        << "  --doubling-family    select every (m,m,2m) signature\n"
        << "  --threads N          evolution/triad worker threads\n"
        << "  --evolve-steps N     evolve input and freeze normalization at input\n"
        << "  --nu X               viscosity for optional evolution\n"
        << "  --dt X               RK4 step for optional evolution\n"
        << "  --backend NAME       direct oracle, fft, or auto (default auto)\n"
        << "  --certificate PATH   write English JSON decomposition\n";
}

int LocalSldSignatureBlockCli::run(
    const LocalSldSignatureBlockCliOptions& options,
    std::ostream& out) {
    SpectralGalerkin galerkin;
    galerkin.configure(options.backend, options.threads);
    const SpectralDynamics dynamics(galerkin);
    SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralReal initial_frequency = 0.0L;
    SpectralReal initial_ep_shift = 0.0L;
    if (options.evolve_steps > 0) {
        const LocalQuarticClosureObjectiveValue initial =
            LocalQuarticClosureObjective(
                dynamics, TriadPartition::local).evaluate(state);
        initial_frequency = initial.initial_frequency;
        initial_ep_shift = initial.initial_ep_shift;
        for (int step = 0; step < options.evolve_steps; ++step) {
            dynamics.rk4_step(
                state, options.viscosity, options.time_step);
        }
    }
    LocalSldSignatureBlockReport report =
        LocalSldSignatureBlock::analyze(
            dynamics, state, options.squared_lengths,
            options.equal_low_doubling_family,
            initial_frequency, initial_ep_shift);
    report.evolved_steps = options.evolve_steps;
    report.evolved_time = options.time_step *
        static_cast<SpectralReal>(options.evolve_steps);
    report.viscosity = options.viscosity;
    write_certificate(report, options);
    out << std::setprecision(12)
        << "local SLD signature block "
        << (report.equal_low_doubling_family
                ? "equal-low doubling family (m,m,2m)"
                : "exact signature")
        << " ("
        << report.squared_lengths[0] << ','
        << report.squared_lengths[1] << ','
        << report.squared_lengths[2] << ")\n"
        << "full_ratio="
        << static_cast<double>(report.full_signed_sld_ratio)
        << " dominant_closed="
        << static_cast<double>(report.dominant_closed_sld_ratio)
        << " remainder_closed="
        << static_cast<double>(report.remainder_closed_sld_ratio)
        << " cross=" << static_cast<double>(report.cross_sld_ratio)
        << " dominant_absolute_fraction="
        << static_cast<double>(report.dominant_absolute_fraction)
        << " reconstruction_error="
        << static_cast<double>(report.quotient_reconstruction_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.finite ? 0 : 2;
}

}  // namespace lemma
