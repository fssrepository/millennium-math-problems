#include "local_sld_cyclic_ansatz.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <cmath>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

SpectralReal pairing(const SpectralIncrement& left,
                     const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        result += std::real(dot_hermitian(left[mode], right[mode]));
    }
    return result;
}

SpectralState cyclic_axis_state() {
    std::mt19937_64 generator(1);
    SpectralState state = SpectralStateFactory::random(2, generator);
    for (ComplexVector& value : state.velocity) {
        value = {};
    }
    auto set_pair = [&](WaveVector wave, std::size_t component) {
        ComplexVector value{};
        value[component] = 1.0L;
        state.velocity[state.index.at(wave)] = value;
        state.velocity[state.index.at(-wave)] = conjugate(value);
    };
    set_pair(WaveVector{1, 0, 0}, 2);
    set_pair(WaveVector{0, 1, 0}, 0);
    set_pair(WaveVector{0, 0, 1}, 1);
    SpectralStateOps::normalize_energy(state);
    return state;
}

SpectralState response_state(const SpectralDynamics& dynamics,
                             const SpectralState& axis) {
    SpectralState response = axis;
    response.velocity = dynamics.advection_direct_partition(
        axis, TriadPartition::local);
    dynamics.enforce_constraints(response);
    SpectralStateOps::normalize_energy(response);
    return response;
}

SpectralState mix(const SpectralState& axis,
                  const SpectralState& response,
                  SpectralReal angle) {
    SpectralState state = axis;
    const SpectralReal left = std::cos(angle);
    const SpectralReal right = std::sin(angle);
    for (std::size_t mode = 0; mode < state.velocity.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            state.velocity[mode][component] =
                left * axis.velocity[mode][component] +
                right * response.velocity[mode][component];
        }
    }
    SpectralStateOps::normalize_energy(state);
    return state;
}

SpectralReal tangent_gradient_norm(
    const LocalQuarticClosureObjective& objective,
    const SpectralState& state) {
    SpectralIncrement gradient =
        objective.signed_local_sld_ratio_gradient(state);
    const SpectralReal energy = SpectralStateOps::energy(state);
    const SpectralReal radial = pairing(gradient, state.velocity) / energy;
    for (std::size_t mode = 0; mode < gradient.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            gradient[mode][component] -=
                radial * state.velocity[mode][component];
        }
    }
    return std::sqrt(std::max(0.0L, pairing(gradient, gradient)));
}

void write_certificate(const LocalSldCyclicAnsatzReport& report,
                       const LocalSldCyclicAnsatzOptions& options) {
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot write cyclic ansatz certificate");
    }
    const auto& value = report.value;
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-cyclic-ansatz-v1\",\n"
        << "  \"ansatz\": \"normalized cyclic axis shear plus its normalized local quadratic advection response\",\n"
        << "  \"angle\": " << static_cast<double>(report.angle) << ",\n"
        << "  \"axis_energy_fraction\": "
        << static_cast<double>(report.axis_energy_fraction) << ",\n"
        << "  \"response_energy_fraction\": "
        << static_cast<double>(report.response_energy_fraction) << ",\n"
        << "  \"basis_inner_product\": "
        << static_cast<double>(report.basis_inner_product) << ",\n"
        << "  \"pure_axis_energy\": "
        << static_cast<double>(report.pure_axis_value.energy) << ",\n"
        << "  \"pure_axis_enstrophy\": "
        << static_cast<double>(report.pure_axis_value.enstrophy) << ",\n"
        << "  \"pure_axis_palinstrophy\": "
        << static_cast<double>(report.pure_axis_value.palinstrophy) << ",\n"
        << "  \"pure_axis_stretching\": "
        << static_cast<double>(report.pure_axis_value.signed_stretching)
        << ",\n"
        << "  \"pure_axis_two_entry_bracket\": "
        << static_cast<double>(
               report.pure_axis_value.signed_two_entry_bracket)
        << ",\n"
        << "  \"pure_axis_signed_constant_ratio\": "
        << static_cast<double>(
               report.pure_axis_value.signed_constant_ratio)
        << ",\n"
        << "  \"pure_axis_absolute_constant_ratio\": "
        << static_cast<double>(report.pure_axis_value.constant_ratio)
        << ",\n"
        << "  \"pure_axis_identity_error\": "
        << static_cast<double>(report.pure_axis_identity_error) << ",\n"
        << "  \"pure_response_signed_constant_ratio\": "
        << static_cast<double>(
               report.pure_response_value.signed_constant_ratio)
        << ",\n"
        << "  \"signed_local_sld_ratio\": "
        << static_cast<double>(value.signed_local_sld_ratio) << ",\n"
        << "  \"closure_constant_ratio\": "
        << static_cast<double>(value.constant_ratio) << ",\n"
        << "  \"normalized_stretching_ratio\": "
        << static_cast<double>(value.normalized_stretching_ratio) << ",\n"
        << "  \"signed_shape_factor\": "
        << static_cast<double>(value.signed_shape_factor) << ",\n"
        << "  \"factorization_relative_error\": "
        << static_cast<double>(value.factorization_relative_error) << ",\n"
        << "  \"signed_stretching\": "
        << static_cast<double>(value.signed_stretching) << ",\n"
        << "  \"signed_two_entry_bracket\": "
        << static_cast<double>(value.signed_two_entry_bracket) << ",\n"
        << "  \"energy\": " << static_cast<double>(value.energy) << ",\n"
        << "  \"enstrophy\": "
        << static_cast<double>(value.enstrophy) << ",\n"
        << "  \"palinstrophy\": "
        << static_cast<double>(value.palinstrophy) << ",\n"
        << "  \"projected_full_gradient_norm\": "
        << static_cast<double>(report.projected_gradient_norm) << ",\n"
        << "  \"coarse_samples\": " << report.coarse_samples << ",\n"
        << "  \"refinement_iterations\": "
        << report.refinement_iterations << ",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"candidate_lemma_proved\": false,\n"
        << "  \"finite_search_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldCyclicAnsatzReport LocalSldCyclicAnsatz::optimize(
    const LocalSldCyclicAnsatzOptions& options) {
    if (options.coarse_samples < 32 ||
        options.refinement_iterations < 1) {
        throw std::invalid_argument("invalid cyclic ansatz resolution");
    }
    SpectralGalerkin galerkin;
    galerkin.configure("direct", 1);
    const SpectralDynamics dynamics(galerkin);
    const LocalQuarticClosureObjective objective(dynamics);
    const SpectralState axis = cyclic_axis_state();
    const SpectralState response = response_state(dynamics, axis);
    const SpectralReal period =
        2.0L * std::numbers::pi_v<SpectralReal>;
    const SpectralReal spacing = period /
        static_cast<SpectralReal>(options.coarse_samples);
    SpectralReal best_angle = 0.0L;
    SpectralReal best_value = -1.0L;
    for (int sample = 0; sample < options.coarse_samples; ++sample) {
        const SpectralReal angle = spacing *
            static_cast<SpectralReal>(sample);
        const SpectralReal value = objective.evaluate(
            mix(axis, response, angle)).signed_local_sld_ratio;
        if (value > best_value) {
            best_value = value;
            best_angle = angle;
        }
    }
    SpectralReal lower = best_angle - spacing;
    SpectralReal upper = best_angle + spacing;
    const SpectralReal golden =
        (std::sqrt(5.0L) - 1.0L) / 2.0L;
    SpectralReal left = upper - golden * (upper - lower);
    SpectralReal right = lower + golden * (upper - lower);
    auto value_at = [&](SpectralReal angle) {
        return objective.evaluate(
            mix(axis, response, angle)).signed_local_sld_ratio;
    };
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
    LocalSldCyclicAnsatzReport report;
    report.angle = 0.5L * (lower + upper);
    report.axis_energy_fraction =
        std::cos(report.angle) * std::cos(report.angle);
    report.response_energy_fraction =
        std::sin(report.angle) * std::sin(report.angle);
    report.basis_inner_product = pairing(
        axis.velocity, response.velocity);
    report.pure_axis_value = objective.evaluate(axis);
    report.pure_response_value = objective.evaluate(response);
    report.pure_axis_identity_error = std::max({
        std::abs(report.pure_axis_value.energy - 1.0L),
        std::abs(report.pure_axis_value.enstrophy - 1.0L),
        std::abs(report.pure_axis_value.palinstrophy - 1.0L),
        std::abs(report.pure_axis_value.signed_stretching),
        std::abs(
            report.pure_axis_value.signed_two_entry_bracket + 1.0L / 3.0L),
        std::abs(
            report.pure_axis_value.signed_constant_ratio + 1.0L / 3.0L)});
    report.state = mix(axis, response, report.angle);
    report.value = objective.evaluate(report.state);
    report.projected_gradient_norm = tangent_gradient_norm(
        objective, report.state);
    report.coarse_samples = options.coarse_samples;
    report.refinement_iterations = options.refinement_iterations;
    return report;
}

LocalSldCyclicAnsatzOptions LocalSldCyclicAnsatzCli::parse(
    int argc, char** argv, int first) {
    LocalSldCyclicAnsatzOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--samples") {
            options.coarse_samples = std::stoi(next(index, name));
        } else if (name == "--refinements") {
            options.refinement_iterations = std::stoi(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--state") {
            options.state_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-sld-ansatz option: " + name);
        }
    }
    if (options.coarse_samples < 32 ||
        options.refinement_iterations < 1 ||
        options.certificate_path.empty() || options.state_path.empty()) {
        throw std::invalid_argument(
            "local-sld-ansatz requires resolution, --certificate, and --state");
    }
    return options;
}

void LocalSldCyclicAnsatzCli::print_help(std::ostream& out) {
    out << "Local SLD cyclic two-basis ansatz options:\n"
        << "  --samples N          coarse periodic angle samples\n"
        << "  --refinements N      golden-section refinement iterations\n"
        << "  --certificate PATH   write English JSON certificate\n"
        << "  --state PATH         write optimized Fourier state\n";
}

int LocalSldCyclicAnsatzCli::run(
    const LocalSldCyclicAnsatzOptions& options,
    std::ostream& out) {
    const LocalSldCyclicAnsatzReport report =
        LocalSldCyclicAnsatz::optimize(options);
    SpectralStateWriter::write_tsv(
        options.state_path, report.state,
        "cyclic shear plus quadratic-response ansatz; candidate_lemma_proved=false");
    write_certificate(report, options);
    out << std::setprecision(12)
        << "cyclic ansatz angle=" << static_cast<double>(report.angle)
        << " axis_energy="
        << static_cast<double>(report.axis_energy_fraction)
        << " response_energy="
        << static_cast<double>(report.response_energy_fraction)
        << " direct_SLD_ratio="
        << static_cast<double>(report.value.signed_local_sld_ratio)
        << " pure_axis_c="
        << static_cast<double>(
               report.pure_axis_value.signed_constant_ratio)
        << " pure_axis_error="
        << static_cast<double>(report.pure_axis_identity_error)
        << " full_gradient_norm="
        << static_cast<double>(report.projected_gradient_norm) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return 0;
}

}  // namespace lemma
