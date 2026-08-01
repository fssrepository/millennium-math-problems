#pragma once

#include "spectral_galerkin.hpp"
#include "triad_partition.hpp"
#include "spectral_state.hpp"

namespace lemma {

class SpectralDynamics {
public:
    explicit SpectralDynamics(const SpectralGalerkin& configuration);

    [[nodiscard]] SpectralIncrement advection_direct(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement advection_direct_partition(
        const SpectralState& state, TriadPartition partition) const;
    [[nodiscard]] SpectralIncrement advection_fft(
        const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement advection_jvp_direct(
        const SpectralState& state,
        const SpectralIncrement& direction) const;
    [[nodiscard]] SpectralIncrement advection_jvp_fft(
        const SpectralState& state,
        const SpectralIncrement& direction) const;
    [[nodiscard]] SpectralIncrement advection_vjp_direct(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent) const;
    [[nodiscard]] SpectralIncrement advection_vjp_direct_partition(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent,
        TriadPartition partition) const;
    [[nodiscard]] SpectralIncrement advection_vjp_fft(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent) const;
    [[nodiscard]] SpectralIncrement advection(const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement advection_jvp(
        const SpectralState& state,
        const SpectralIncrement& direction) const;
    [[nodiscard]] SpectralIncrement advection_vjp(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent) const;
    [[nodiscard]] SpectralIncrement rhs(const SpectralState& state,
                                        SpectralReal viscosity) const;
    [[nodiscard]] SpectralIncrement rhs_jvp(
        const SpectralState& state, const SpectralIncrement& direction,
        SpectralReal viscosity) const;
    [[nodiscard]] SpectralIncrement rhs_vjp(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent,
        SpectralReal viscosity) const;
    [[nodiscard]] SpectralState add_increment(
        const SpectralState& base, const SpectralIncrement& increment,
        SpectralReal factor) const;
    void enforce_constraints(SpectralState& state) const;
    void rk4_step(SpectralState& state, SpectralReal viscosity,
                  SpectralReal time_step) const;
    [[nodiscard]] SpectralIncrement rk4_jvp(
        const SpectralState& state, const SpectralIncrement& direction,
        SpectralReal viscosity, SpectralReal time_step) const;
    [[nodiscard]] SpectralIncrement rk4_vjp(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent,
        SpectralReal viscosity, SpectralReal time_step) const;

private:
    void project_increment(SpectralIncrement& increment,
                           const SpectralState& layout) const;
    const SpectralGalerkin& configuration_;
};

}  // namespace lemma
