#include "triad_verifier.hpp"

#include "triad_partition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lemma {
namespace {

using TriadKey = std::array<WaveVector, 3>;

}  // namespace

struct TriadVerifier::InteractionAnalysis {
    std::vector<ComplexVector> advection;
    std::vector<SpectralReal> local_energy_transfer;
    std::vector<SpectralReal> nonlocal_energy_transfer;
    SpectralReal local_absolute_transfer = 0.0L;
    SpectralReal nonlocal_absolute_transfer = 0.0L;
    SpectralReal maximum_detailed_triad_residual = 0.0L;
    SpectralReal maximum_relative_detailed_triad_residual = 0.0L;
};

struct TriadVerifier::TriadMeasurement {
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal advection_norm = 0.0L;
    SpectralReal energy_pairing = 0.0L;
    SpectralReal vortex_stretching = 0.0L;
    SpectralReal divergence_residual = 0.0L;
    SpectralReal reality_residual = 0.0L;
    SpectralReal classical_ratio = 0.0L;
    SpectralReal detailed_triad_residual = 0.0L;
    SpectralReal relative_detailed_triad_residual = 0.0L;
    SpectralReal local_absolute_transfer = 0.0L;
    SpectralReal nonlocal_absolute_transfer = 0.0L;
    SpectralReal maximum_cumulative_flux = 0.0L;
    SpectralReal maximum_local_cumulative_flux = 0.0L;
    SpectralReal maximum_nonlocal_cumulative_flux = 0.0L;
    SpectralReal flux_partition_residual = 0.0L;
};

TriadVerifier::InteractionAnalysis TriadVerifier::analyze_interactions(
    const SpectralState& state) {
    InteractionAnalysis result;
    result.advection.resize(state.waves.size());
    result.local_energy_transfer.assign(state.waves.size(), 0.0L);
    result.nonlocal_energy_transfer.assign(state.waves.size(), 0.0L);
    std::map<TriadKey, std::pair<SpectralReal, SpectralReal>> detailed_triads;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, target_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const ComplexVector& up = state.velocity[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[target_index];
        const SpectralComplex coefficient =
            imaginary_unit * wave_dot(q, up);
        const ComplexVector& uq = state.velocity[q_index];
        ComplexVector& value = result.advection[target_index];
        ComplexVector pair_contribution{};
        for (std::size_t direction = 0; direction < 3; ++direction) {
            pair_contribution[direction] = coefficient * uq[direction];
            value[direction] += pair_contribution[direction];
        }

        const SpectralReal transfer = std::real(dot_hermitian(
            state.velocity[target_index], pair_contribution));
        if (TriadPartitioner::is_local(k, p, q)) {
            result.local_energy_transfer[target_index] += transfer;
            result.local_absolute_transfer += std::abs(transfer);
        } else {
            result.nonlocal_energy_transfer[target_index] += transfer;
            result.nonlocal_absolute_transfer += std::abs(transfer);
        }

        TriadKey key{-k, p, q};
        std::sort(key.begin(), key.end());
        auto& [sum, absolute_sum] = detailed_triads[key];
        sum += transfer;
        absolute_sum += std::abs(transfer);
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        result.advection[index] = project_divergence_free(
            state.waves[index], result.advection[index]);
    }

    const SpectralReal total_absolute_transfer =
        result.local_absolute_transfer + result.nonlocal_absolute_transfer;
    for (const auto& [key, cancellation] : detailed_triads) {
        static_cast<void>(key);
        const auto [sum, absolute_sum] = cancellation;
        result.maximum_detailed_triad_residual = std::max(
            result.maximum_detailed_triad_residual, std::abs(sum));
        if (absolute_sum > 1e-14L * total_absolute_transfer) {
            result.maximum_relative_detailed_triad_residual = std::max(
                result.maximum_relative_detailed_triad_residual,
                std::abs(sum) / absolute_sum);
        }
    }
    if (total_absolute_transfer > 0.0L) {
        result.maximum_detailed_triad_residual /= total_absolute_transfer;
    }
    return result;
}

TriadVerifier::TriadMeasurement TriadVerifier::measure(
    const SpectralState& state) {
    const InteractionAnalysis interactions = analyze_interactions(state);
    TriadMeasurement result;
    result.detailed_triad_residual =
        interactions.maximum_detailed_triad_residual;
    result.relative_detailed_triad_residual =
        interactions.maximum_relative_detailed_triad_residual;
    result.local_absolute_transfer = interactions.local_absolute_transfer;
    result.nonlocal_absolute_transfer = interactions.nonlocal_absolute_transfer;
    int maximum_radius = 0;
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(wave));
        maximum_radius = std::max(
            maximum_radius,
            static_cast<int>(std::ceil(std::sqrt(wave2))));
        const SpectralReal velocity2 = std::real(dot_hermitian(
            state.velocity[index], state.velocity[index]));
        const SpectralReal nonlinear2 = std::real(dot_hermitian(
            interactions.advection[index], interactions.advection[index]));
        const SpectralReal pairing = std::real(dot_hermitian(
            state.velocity[index], interactions.advection[index]));
        result.energy += velocity2;
        result.enstrophy += wave2 * velocity2;
        result.palinstrophy += wave2 * wave2 * velocity2;
        result.advection_norm += nonlinear2;
        result.energy_pairing += pairing;
        result.vortex_stretching += wave2 * pairing;
        result.divergence_residual = std::max(
            result.divergence_residual,
            std::abs(wave_dot(wave, state.velocity[index])));
        const auto opposite = state.index.find(-wave);
        if (opposite != state.index.end()) {
            for (std::size_t direction = 0; direction < 3; ++direction) {
                result.reality_residual = std::max(
                    result.reality_residual,
                    std::abs(
                        state.velocity[opposite->second][direction] -
                        std::conj(state.velocity[index][direction])));
            }
        }
    }
    const SpectralReal denominator =
        std::pow(result.enstrophy, 0.75L) *
        std::pow(result.palinstrophy, 0.75L);
    if (denominator > 0.0L) {
        result.classical_ratio =
            std::abs(result.vortex_stretching) / denominator;
    }

    for (int cutoff = 1; cutoff < maximum_radius; ++cutoff) {
        SpectralReal local_flux = 0.0L;
        SpectralReal nonlocal_flux = 0.0L;
        const SpectralInteger cutoff2 =
            static_cast<SpectralInteger>(cutoff) * cutoff;
        for (std::size_t index = 0; index < state.waves.size(); ++index) {
            if (norm_squared(state.waves[index]) > cutoff2) {
                local_flux -= interactions.local_energy_transfer[index];
                nonlocal_flux -= interactions.nonlocal_energy_transfer[index];
            }
        }
        const SpectralReal total_flux = local_flux + nonlocal_flux;
        result.maximum_cumulative_flux = std::max(
            result.maximum_cumulative_flux, std::abs(total_flux));
        result.maximum_local_cumulative_flux = std::max(
            result.maximum_local_cumulative_flux, std::abs(local_flux));
        result.maximum_nonlocal_cumulative_flux = std::max(
            result.maximum_nonlocal_cumulative_flux,
            std::abs(nonlocal_flux));
        result.flux_partition_residual = std::max(
            result.flux_partition_residual,
            std::abs(total_flux - local_flux - nonlocal_flux));
    }
    return result;
}

TriadCertificate TriadVerifier::analyze(
    int cutoff, int samples, std::uint64_t seed) {
    if (samples < 1 || samples > 100000) {
        throw std::invalid_argument(
            "--triad-samples must be between 1 and 100000");
    }
    std::mt19937_64 generator(seed);
    TriadCertificate certificate;
    certificate.samples = samples;
    for (int sample = 0; sample < samples; ++sample) {
        const SpectralState state =
            SpectralStateFactory::random(cutoff, generator);
        certificate.modes = static_cast<int>(state.waves.size());
        const TriadMeasurement measurement = measure(state);
        const SpectralReal scale = std::sqrt(std::max(
            0.0L, measurement.energy * measurement.advection_norm));
        const SpectralReal normalized_energy_residual =
            scale > 0.0L
                ? std::abs(measurement.energy_pairing) / scale
                : 0.0L;
        certificate.maximum_normalized_energy_residual = std::max(
            certificate.maximum_normalized_energy_residual,
            normalized_energy_residual);
        certificate.maximum_divergence_residual = std::max(
            certificate.maximum_divergence_residual,
            measurement.divergence_residual);
        certificate.maximum_reality_residual = std::max(
            certificate.maximum_reality_residual,
            measurement.reality_residual);
        certificate.maximum_classical_ratio = std::max(
            certificate.maximum_classical_ratio,
            measurement.classical_ratio);
        certificate.maximum_vortex_stretching = std::max(
            certificate.maximum_vortex_stretching,
            std::abs(measurement.vortex_stretching));
        certificate.maximum_detailed_triad_residual = std::max(
            certificate.maximum_detailed_triad_residual,
            measurement.detailed_triad_residual);
        certificate.maximum_relative_detailed_triad_residual = std::max(
            certificate.maximum_relative_detailed_triad_residual,
            measurement.relative_detailed_triad_residual);
        const SpectralReal absolute_transfer =
            measurement.local_absolute_transfer +
            measurement.nonlocal_absolute_transfer;
        if (absolute_transfer > 0.0L) {
            certificate.maximum_nonlocal_absolute_fraction = std::max(
                certificate.maximum_nonlocal_absolute_fraction,
                measurement.nonlocal_absolute_transfer / absolute_transfer);
            certificate.maximum_flux_efficiency = std::max(
                certificate.maximum_flux_efficiency,
                measurement.maximum_cumulative_flux / absolute_transfer);
        }
        certificate.maximum_local_cumulative_flux = std::max(
            certificate.maximum_local_cumulative_flux,
            measurement.maximum_local_cumulative_flux);
        certificate.maximum_nonlocal_cumulative_flux = std::max(
            certificate.maximum_nonlocal_cumulative_flux,
            measurement.maximum_nonlocal_cumulative_flux);
        certificate.maximum_flux_partition_residual = std::max(
            certificate.maximum_flux_partition_residual,
            measurement.flux_partition_residual);
        certificate.nonzero_vortex_stretching_seen =
            certificate.nonzero_vortex_stretching_seen ||
            std::abs(measurement.vortex_stretching) > 1e-16L;
    }
    return certificate;
}

}  // namespace lemma
