#include "local_sld_remainder_double_square.hpp"

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

SpectralReal pairing(const SpectralIncrement& left,
                     const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "remainder double-square layout mismatch");
    }
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralIncrement laplacian_weight(
    const SpectralState& state,
    const SpectralIncrement& source,
    SpectralReal exponent) {
    if (source.size() != state.waves.size()) {
        throw std::invalid_argument(
            "remainder double-square Laplacian layout mismatch");
    }
    SpectralIncrement result = source;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal weight = std::pow(wave2, exponent);
        for (SpectralComplex& component : result[mode]) {
            component *= weight;
        }
    }
    return result;
}

SpectralIncrement linear_combination(
    const SpectralIncrement& left,
    SpectralReal left_scale,
    const SpectralIncrement& right,
    SpectralReal right_scale) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "remainder double-square combination mismatch");
    }
    SpectralIncrement result = left;
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                left_scale * left[mode][component] +
                right_scale * right[mode][component];
        }
    }
    return result;
}

SpectralReal relative_error(SpectralReal computed,
                            SpectralReal expected) {
    return std::abs(computed - expected) /
        std::max({std::abs(computed), std::abs(expected), 1e-30L});
}

SpectralReal increment_relative_error(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    if (left.size() != right.size()) {
        throw std::invalid_argument(
            "remainder double-square residual layout mismatch");
    }
    SpectralReal residual2 = 0.0L;
    SpectralReal scale2 = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            residual2 += std::norm(
                left[mode][component] - right[mode][component]);
            scale2 += std::norm(left[mode][component]);
            scale2 += std::norm(right[mode][component]);
        }
    }
    return std::sqrt(residual2 / std::max(scale2, 1e-30L));
}

void write_json(
    const LocalSldRemainderDoubleSquareReport& report,
    const LocalSldRemainderDoubleSquareCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-remainder-double-square-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"selection\": \"local signatures excluding (m,m,2m)\",\n"
        << "  \"identity\": \"K+G = -||A^(1/2)(B-cAu-(1/2)A^(-1)D)||_2^2 + S^2/(2Z) + c^2||A^(3/2)u||_2^2 + c<Au,D> + (1/4)||A^(-1/2)D||_2^2\",\n"
        << "  \"commutator\": \"D=B(u,Au)-[x -> B(x,u)]^*Au\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"enstrophy\": "
        << static_cast<double>(report.enstrophy) << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(report.palinstrophy) << ",\n"
        << "  \"stretching\": "
        << static_cast<double>(report.stretching) << ",\n"
        << "  \"palinstrophy_cross\": "
        << static_cast<double>(report.palinstrophy_cross) << ",\n"
        << "  \"projection_coefficient\": "
        << static_cast<double>(report.projection_coefficient) << ",\n"
        << "  \"full_bracket\": "
        << static_cast<double>(report.full_bracket) << ",\n"
        << "  \"target_scale\": "
        << static_cast<double>(report.target_scale) << ",\n"
        << "  \"signed_target_ratio\": "
        << static_cast<double>(report.signed_target_ratio) << ",\n"
        << "  \"first_square_norm2\": "
        << static_cast<double>(report.first_square_norm2) << ",\n"
        << "  \"enstrophy_normalization\": "
        << static_cast<double>(
               report.enstrophy_normalization) << ",\n"
        << "  \"projected_h3_correction\": "
        << static_cast<double>(
               report.projected_h3_correction) << ",\n"
        << "  \"commutator_pairing\": "
        << static_cast<double>(report.commutator_pairing) << ",\n"
        << "  \"projected_commutator_pairing\": "
        << static_cast<double>(
               report.projected_commutator_pairing) << ",\n"
        << "  \"commutator_hminus1_norm2\": "
        << static_cast<double>(
               report.commutator_hminus1_norm2) << ",\n"
        << "  \"second_square_norm2\": "
        << static_cast<double>(report.second_square_norm2) << ",\n"
        << "  \"completed_upper_envelope\": "
        << static_cast<double>(
               report.completed_upper_envelope) << ",\n"
        << "  \"upper_envelope_target_ratio\": "
        << static_cast<double>(
               report.upper_envelope_target_ratio) << ",\n"
        << "  \"negative_second_square_target_ratio\": "
        << static_cast<double>(
               report.negative_second_square_target_ratio) << ",\n"
        << "  \"first_completion_error\": "
        << static_cast<double>(report.first_completion_error) << ",\n"
        << "  \"commutator_reconstruction_error\": "
        << static_cast<double>(
               report.commutator_reconstruction_error) << ",\n"
        << "  \"stretching_vjp_identity\": \"D=-[dB(u,u)]^*Au and grad(S)=AB-D\",\n"
        << "  \"stretching_vjp_reconstruction_error\": "
        << static_cast<double>(
               report.stretching_vjp_reconstruction_error) << ",\n"
        << "  \"second_completion_error\": "
        << static_cast<double>(report.second_completion_error) << ",\n"
        << "  \"exact_identity\": "
        << (report.exact_identity ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_upper_bound_proved\": false,\n"
        << "  \"remaining_requirement\": \"bound the three nonnegative-envelope terms at Z^(5/4)P^(3/4), exploiting the commutator D rather than unsigned dense incidence\",\n"
        << "  \"finite_identity_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldRemainderDoubleSquareReport
LocalSldRemainderDoubleSquare::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state) {
    const TriadSelection selection =
        TriadSelection::local_without_equal_low_doubling();
    const LocalQuarticClosureObjectiveValue objective =
        LocalQuarticClosureObjective(dynamics, selection).evaluate(state);
    const SpectralIncrement au = laplacian_weight(
        state, state.velocity, 1.0L);
    const SpectralIncrement b =
        dynamics.advection_direct_partition(state, selection);
    const SpectralIncrement transported_au =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, au, selection);
    const BilinearAdvectionCotangents nested_cotangents =
        dynamics.advection_bilinear_vjp_direct_partition(
            state, state.velocity, state.velocity, au, selection);
    const SpectralIncrement d = linear_combination(
        transported_au, 1.0L,
        nested_cotangents.advecting, -1.0L);
    SpectralIncrement negative_stretching_vjp =
        dynamics.advection_vjp_direct_partition(
            state, au, selection);
    for (ComplexVector& mode : negative_stretching_vjp) {
        for (SpectralComplex& component : mode) {
            component = -component;
        }
    }
    const SpectralIncrement inverse_d = laplacian_weight(
        state, d, -1.0L);

    LocalSldRemainderDoubleSquareReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.enstrophy = objective.enstrophy;
    report.palinstrophy = objective.palinstrophy;
    report.stretching = objective.signed_stretching;
    report.palinstrophy_cross = objective.palinstrophy_cross;
    report.projection_coefficient = 3.0L * report.stretching /
        (4.0L * report.palinstrophy);
    report.full_bracket = objective.signed_two_entry_bracket;
    report.target_scale = objective.lqc3_target_scale;
    report.signed_target_ratio = objective.signed_lqc3_target_ratio;

    const SpectralIncrement w = linear_combination(
        b, 1.0L, au, -report.projection_coefficient);
    const SpectralIncrement aw = laplacian_weight(state, w, 1.0L);
    report.first_square_norm2 = pairing(w, aw);
    report.enstrophy_normalization = report.stretching *
        report.stretching / (2.0L * report.enstrophy);
    const SpectralReal h3_squared = pairing(
        au, laplacian_weight(state, au, 1.0L));
    report.projected_h3_correction =
        report.projection_coefficient *
        report.projection_coefficient * h3_squared;
    report.commutator_pairing = pairing(b, d);
    report.projected_commutator_pairing =
        report.projection_coefficient * pairing(au, d);
    report.commutator_hminus1_norm2 = pairing(d, inverse_d);

    const SpectralIncrement twice_completed = linear_combination(
        w, 1.0L, inverse_d, -0.5L);
    report.second_square_norm2 = pairing(
        twice_completed,
        laplacian_weight(state, twice_completed, 1.0L));
    report.completed_upper_envelope =
        report.enstrophy_normalization +
        report.projected_h3_correction +
        report.projected_commutator_pairing +
        0.25L * report.commutator_hminus1_norm2;
    if (report.target_scale > 0.0L) {
        report.upper_envelope_target_ratio =
            report.completed_upper_envelope / report.target_scale;
        report.negative_second_square_target_ratio =
            -report.second_square_norm2 / report.target_scale;
    }

    const SpectralReal first_reconstruction =
        -report.first_square_norm2 +
        report.enstrophy_normalization +
        report.projected_h3_correction +
        report.commutator_pairing;
    const SpectralReal direct_nonprojected =
        objective.negative_commutator_pairing +
        objective.advecting_slot + pairing(
            b, laplacian_weight(state, b, 1.0L));
    report.first_completion_error = relative_error(
        first_reconstruction, report.full_bracket);
    report.commutator_reconstruction_error = relative_error(
        report.commutator_pairing, direct_nonprojected);
    report.stretching_vjp_reconstruction_error =
        increment_relative_error(d, negative_stretching_vjp);
    report.second_completion_error = relative_error(
        -report.second_square_norm2 +
            report.completed_upper_envelope,
        report.full_bracket);
    report.exact_identity =
        report.first_completion_error < 1e-13L &&
        report.commutator_reconstruction_error < 1e-13L &&
        report.stretching_vjp_reconstruction_error < 1e-13L &&
        report.second_completion_error < 1e-13L &&
        report.commutator_hminus1_norm2 >= -1e-14L &&
        report.first_square_norm2 >= -1e-14L &&
        report.second_square_norm2 >= -1e-14L;
    return report;
}

LocalSldRemainderDoubleSquareCliOptions
LocalSldRemainderDoubleSquareCli::parse(
    int argc, char** argv, int first) {
    LocalSldRemainderDoubleSquareCliOptions options;
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
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown remainder-double-square option: " + name);
        }
    }
    if (options.state_path.empty() ||
        options.certificate_path.empty() || options.threads < 1) {
        throw std::invalid_argument(
            "remainder-double-square requires state, certificate and positive threads");
    }
    return options;
}

void LocalSldRemainderDoubleSquareCli::print_help(std::ostream& out) {
    out << "Local SLD remainder double-square options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON identity\n"
        << "  --threads N           direct-kernel workers\n";
}

int LocalSldRemainderDoubleSquareCli::run(
    const LocalSldRemainderDoubleSquareCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldRemainderDoubleSquareReport report =
        LocalSldRemainderDoubleSquare::analyze(dynamics, state);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write remainder double-square certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12)
        << "remainder double square cutoff=" << report.cutoff
        << " signed_lqc3="
        << static_cast<double>(report.signed_target_ratio)
        << " upper_envelope="
        << static_cast<double>(report.upper_envelope_target_ratio)
        << " negative_square="
        << static_cast<double>(report.negative_second_square_target_ratio)
        << " identity_error="
        << static_cast<double>(report.second_completion_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_identity ? 0 : 2;
}

}  // namespace lemma
