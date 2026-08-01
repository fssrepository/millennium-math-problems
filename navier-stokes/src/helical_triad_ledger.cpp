#include "helical_triad_ledger.hpp"

#include "triad_ledger.hpp"
#include "triad_partition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace lemma {
namespace {

using RealVector = std::array<SpectralReal, 3>;

RealVector cross(const RealVector& left, const RealVector& right) {
    return {left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0]};
}

SpectralReal norm(const RealVector& value) {
    return std::sqrt(value[0] * value[0] + value[1] * value[1] +
                     value[2] * value[2]);
}

std::array<ComplexVector, 2> helical_components(
    const WaveVector& wave, const ComplexVector& velocity) {
    const SpectralReal wave_norm = std::sqrt(
        static_cast<SpectralReal>(norm_squared(wave)));
    const RealVector direction{
        static_cast<SpectralReal>(wave.x) / wave_norm,
        static_cast<SpectralReal>(wave.y) / wave_norm,
        static_cast<SpectralReal>(wave.z) / wave_norm};
    std::size_t reference_index = 0;
    if (std::abs(direction[1]) < std::abs(direction[reference_index])) {
        reference_index = 1;
    }
    if (std::abs(direction[2]) < std::abs(direction[reference_index])) {
        reference_index = 2;
    }
    RealVector reference{};
    reference[reference_index] = 1.0L;
    RealVector first = cross(direction, reference);
    const SpectralReal first_norm = norm(first);
    for (SpectralReal& component : first) {
        component /= first_norm;
    }
    const RealVector second = cross(direction, first);
    const SpectralReal inverse_sqrt_two = 1.0L / std::sqrt(2.0L);
    std::array<ComplexVector, 2> basis{};
    for (std::size_t sign_index = 0; sign_index < 2; ++sign_index) {
        const SpectralReal sign = sign_index == 0 ? -1.0L : 1.0L;
        for (std::size_t component = 0; component < 3; ++component) {
            basis[sign_index][component] = inverse_sqrt_two *
                SpectralComplex{first[component], sign * second[component]};
        }
    }
    std::array<ComplexVector, 2> result{};
    for (std::size_t sign_index = 0; sign_index < 2; ++sign_index) {
        const SpectralComplex coefficient = dot_hermitian(
            basis[sign_index], velocity);
        for (std::size_t component = 0; component < 3; ++component) {
            result[sign_index][component] =
                coefficient * basis[sign_index][component];
        }
    }
    return result;
}

SpectralReal interaction_stretching(
    const WaveVector& advected_wave, SpectralInteger target_wave2,
    const ComplexVector& advecting, const ComplexVector& advected,
    const ComplexVector& target) {
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    const SpectralComplex coefficient = imaginary_unit *
        wave_dot(advected_wave, advecting);
    ComplexVector contribution{};
    for (std::size_t component = 0; component < 3; ++component) {
        contribution[component] = coefficient * advected[component];
    }
    return static_cast<SpectralReal>(target_wave2) *
        std::real(dot_hermitian(target, contribution));
}

SpectralReal relative_residual(
    SpectralReal difference, SpectralReal reference) {
    if (reference == 0.0L) {
        return std::abs(difference);
    }
    return std::abs(difference) / std::abs(reference);
}

}  // namespace

SpectralState HelicalTriadLedger::project_helicity(
    const SpectralState& state, int sign) {
    if (sign != -1 && sign != 1) {
        throw std::invalid_argument("helical sign must be -1 or +1");
    }
    SpectralState projected = state;
    const std::size_t sign_index = sign < 0 ? 0U : 1U;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        projected.velocity[mode] = helical_components(
            state.waves[mode], state.velocity[mode])[sign_index];
    }
    return projected;
}

HelicalTriadReport HelicalTriadLedger::analyze(
    const SpectralState& state) {
    HelicalTriadReport report;
    std::array<SpectralIncrement, 2> components{
        SpectralIncrement(state.waves.size()),
        SpectralIncrement(state.waves.size())};
    SpectralReal velocity_norm2 = 0.0L;
    SpectralReal reconstruction_error2 = 0.0L;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const auto split = helical_components(
            state.waves[mode], state.velocity[mode]);
        components[0][mode] = split[0];
        components[1][mode] = split[1];
        const SpectralReal negative_energy = std::real(
            dot_hermitian(split[0], split[0]));
        const SpectralReal positive_energy = std::real(
            dot_hermitian(split[1], split[1]));
        report.negative_helical_energy += negative_energy;
        report.positive_helical_energy += positive_energy;
        report.helicity += std::sqrt(static_cast<SpectralReal>(
            norm_squared(state.waves[mode]))) *
            (positive_energy - negative_energy);
        ComplexVector reconstructed{};
        for (std::size_t component = 0; component < 3; ++component) {
            reconstructed[component] =
                split[0][component] + split[1][component];
            reconstruction_error2 += std::norm(
                reconstructed[component] - state.velocity[mode][component]);
        }
        velocity_norm2 += std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
    }
    report.relative_velocity_reconstruction_residual = std::sqrt(
        reconstruction_error2 /
        std::max(std::numeric_limits<SpectralReal>::min(), velocity_norm2));

    for (std::size_t sector = 0; sector < helical_sector_count; ++sector) {
        report.sectors[sector].advecting_sign =
            (sector & 4U) == 0U ? -1 : 1;
        report.sectors[sector].advected_sign =
            (sector & 2U) == 0U ? -1 : 1;
        report.sectors[sector].target_sign =
            (sector & 1U) == 0U ? -1 : 1;
    }

    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[k_index];
        const bool local = TriadPartitioner::is_local(p, q, k);
        for (std::size_t p_sign = 0; p_sign < 2; ++p_sign) {
            for (std::size_t q_sign = 0; q_sign < 2; ++q_sign) {
                for (std::size_t k_sign = 0; k_sign < 2; ++k_sign) {
                    const std::size_t sector =
                        (p_sign << 2U) | (q_sign << 1U) | k_sign;
                    const SpectralReal stretching = interaction_stretching(
                        q, norm_squared(k), components[p_sign][p_index],
                        components[q_sign][q_index],
                        components[k_sign][k_index]);
                    HelicalSectorRow& row = report.sectors[sector];
                    row.signed_total_stretching += stretching;
                    row.absolute_total_stretching += std::abs(stretching);
                    if (local) {
                        row.signed_local_stretching += stretching;
                        row.absolute_local_stretching += std::abs(stretching);
                    }
                }
            }
        }
    }
    for (std::size_t sector = 0; sector < helical_sector_count; ++sector) {
        const HelicalSectorRow& row = report.sectors[sector];
        report.signed_total_stretching += row.signed_total_stretching;
        report.signed_local_stretching += row.signed_local_stretching;
        const bool homochiral = sector == 0U || sector == 7U;
        if (homochiral) {
            report.homochiral_local_stretching +=
                row.signed_local_stretching;
            report.homochiral_absolute_local_stretching +=
                row.absolute_local_stretching;
        } else {
            report.heterochiral_local_stretching +=
                row.signed_local_stretching;
            report.heterochiral_absolute_local_stretching +=
                row.absolute_local_stretching;
        }
    }
    const TriadLedgerReport ledger = TriadLedger::analyze(state);
    report.relative_total_reconstruction_residual = relative_residual(
        report.signed_total_stretching - ledger.signed_total,
        ledger.signed_total);
    report.relative_local_reconstruction_residual = relative_residual(
        report.signed_local_stretching - ledger.signed_local,
        ledger.signed_local);
    return report;
}

HelicalTriadCertificate HelicalTriadLedger::verify_random(
    int cutoff, int samples, std::uint64_t seed) {
    if (cutoff < 1 || cutoff > 16 || samples < 1 || samples > 100000) {
        throw std::invalid_argument("invalid helical certificate parameters");
    }
    HelicalTriadCertificate certificate;
    certificate.cutoff = cutoff;
    certificate.samples = samples;
    certificate.seed = seed;
    std::mt19937_64 generator(seed);
    constexpr SpectralReal tolerance = 1e-15L;
    for (int sample = 0; sample < samples; ++sample) {
        SpectralState state = SpectralStateFactory::random(cutoff, generator);
        SpectralStateOps::normalize_energy(state);
        const HelicalTriadReport mixed = analyze(state);
        certificate.maximum_velocity_reconstruction_residual = std::max(
            certificate.maximum_velocity_reconstruction_residual,
            mixed.relative_velocity_reconstruction_residual);
        certificate.maximum_total_reconstruction_residual = std::max(
            certificate.maximum_total_reconstruction_residual,
            mixed.relative_total_reconstruction_residual);
        certificate.maximum_local_reconstruction_residual = std::max(
            certificate.maximum_local_reconstruction_residual,
            mixed.relative_local_reconstruction_residual);
        for (const int sign : {-1, 1}) {
            SpectralState pure = project_helicity(state, sign);
            SpectralStateOps::normalize_energy(pure);
            const HelicalTriadReport report = analyze(pure);
            certificate.maximum_velocity_reconstruction_residual = std::max(
                certificate.maximum_velocity_reconstruction_residual,
                report.relative_velocity_reconstruction_residual);
            certificate.maximum_total_reconstruction_residual = std::max(
                certificate.maximum_total_reconstruction_residual,
                report.relative_total_reconstruction_residual);
            certificate.maximum_local_reconstruction_residual = std::max(
                certificate.maximum_local_reconstruction_residual,
                report.relative_local_reconstruction_residual);
            certificate.maximum_pure_heterochiral_absolute_local = std::max(
                certificate.maximum_pure_heterochiral_absolute_local,
                report.heterochiral_absolute_local_stretching);
            certificate.maximum_pure_homochiral_local_stretching = std::max(
                certificate.maximum_pure_homochiral_local_stretching,
                std::abs(report.homochiral_local_stretching));
            certificate.nonzero_pure_homochiral_local_seen =
                certificate.nonzero_pure_homochiral_local_seen ||
                std::abs(report.homochiral_local_stretching) > 1e-12L;
        }
    }
    certificate.all_reconstruction_checks_hold =
        certificate.maximum_velocity_reconstruction_residual < tolerance &&
        certificate.maximum_total_reconstruction_residual < tolerance &&
        certificate.maximum_local_reconstruction_residual < tolerance &&
        certificate.maximum_pure_heterochiral_absolute_local < tolerance;
    return certificate;
}

}  // namespace lemma
