#include "spectral_dynamics.hpp"

#include "spectral_fft_operator.hpp"

#include <cstddef>
#include <stdexcept>

namespace lemma {
namespace {

void require_matching_increment(const SpectralState& state,
                                const SpectralIncrement& increment) {
    if (increment.size() != state.velocity.size()) {
        throw std::invalid_argument(
            "spectral increment does not match state layout");
    }
}

SpectralIncrement scaled_increment(const SpectralIncrement& source,
                                   SpectralReal factor) {
    SpectralIncrement result = source;
    for (ComplexVector& value : result) {
        for (SpectralComplex& component : value) {
            component *= factor;
        }
    }
    return result;
}

void add_scaled_increment(SpectralIncrement& target,
                          const SpectralIncrement& source,
                          SpectralReal factor) {
    if (target.size() != source.size()) {
        throw std::invalid_argument("cannot add mismatched spectral increments");
    }
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += factor * source[mode][component];
        }
    }
}

}  // namespace

SpectralDynamics::SpectralDynamics(const SpectralGalerkin& configuration)
    : configuration_(configuration) {}

SpectralIncrement SpectralDynamics::advection_direct(
    const SpectralState& state) const {
    return advection_direct_partition(state, TriadPartition::all);
}

SpectralIncrement SpectralDynamics::advection_direct_partition(
    const SpectralState& state, TriadPartition partition) const {
    SpectralIncrement advection_result(state.waves.size());
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (!TriadPartitioner::includes(state, interaction, partition)) {
            continue;
        }
        const auto [p_index, q_index, target_index] = interaction;
        const ComplexVector& up = state.velocity[p_index];
        const SpectralComplex coefficient =
            imaginary_unit * wave_dot(state.waves[q_index], up);
        for (std::size_t direction = 0; direction < 3; ++direction) {
            advection_result[target_index][direction] +=
                coefficient * state.velocity[q_index][direction];
        }
    }
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        advection_result[index] = project_divergence_free(
            state.waves[index], advection_result[index]);
    }
    return advection_result;
}

SpectralIncrement SpectralDynamics::advection_fft(
    const SpectralState& state) const {
    return SpectralFftOperator::advection(
        state, configuration_.compute_threads());
}

SpectralIncrement SpectralDynamics::advection(
    const SpectralState& state) const {
    return configuration_.uses_fft(SpectralStateOps::cutoff(state))
               ? advection_fft(state)
               : advection_direct(state);
}

SpectralIncrement SpectralDynamics::advection_jvp(
    const SpectralState& state,
    const SpectralIncrement& direction) const {
    return configuration_.uses_fft(SpectralStateOps::cutoff(state))
               ? advection_jvp_fft(state, direction)
               : advection_jvp_direct(state, direction);
}

SpectralIncrement SpectralDynamics::advection_jvp_direct(
    const SpectralState& state,
    const SpectralIncrement& direction) const {
    require_matching_increment(state, direction);
    SpectralIncrement result(state.waves.size());
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, target_index] = interaction;
        const SpectralComplex base_coefficient =
            imaginary_unit *
            wave_dot(state.waves[q_index], state.velocity[p_index]);
        const SpectralComplex tangent_coefficient =
            imaginary_unit *
            wave_dot(state.waves[q_index], direction[p_index]);
        for (std::size_t component = 0; component < 3; ++component) {
            result[target_index][component] +=
                tangent_coefficient * state.velocity[q_index][component] +
                base_coefficient * direction[q_index][component];
        }
    }
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        result[mode] =
            project_divergence_free(state.waves[mode], result[mode]);
    }
    return result;
}

SpectralIncrement SpectralDynamics::advection_vjp(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent) const {
    return configuration_.uses_fft(SpectralStateOps::cutoff(state))
               ? advection_vjp_fft(state, output_cotangent)
               : advection_vjp_direct(state, output_cotangent);
}

SpectralIncrement SpectralDynamics::advection_vjp_direct(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent) const {
    return advection_vjp_direct_partition(
        state, output_cotangent, TriadPartition::all);
}

SpectralIncrement SpectralDynamics::advection_vjp_direct_partition(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent,
    TriadPartition partition) const {
    require_matching_increment(state, output_cotangent);
    SpectralIncrement cotangent = output_cotangent;
    project_increment(cotangent, state);
    SpectralIncrement result(state.waves.size());
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (!TriadPartitioner::includes(state, interaction, partition)) {
            continue;
        }
        const auto [p_index, q_index, target_index] = interaction;
        const ComplexVector& target_cotangent = cotangent[target_index];
        const SpectralComplex first_coefficient =
            minus_imaginary_unit *
            dot_hermitian(state.velocity[q_index], target_cotangent);
        for (std::size_t component = 0; component < 3; ++component) {
            const SpectralReal wave_component = static_cast<SpectralReal>(
                component == 0   ? state.waves[q_index].x
                : component == 1 ? state.waves[q_index].y
                                 : state.waves[q_index].z);
            result[p_index][component] +=
                wave_component * first_coefficient;
        }
        const SpectralComplex second_coefficient =
            minus_imaginary_unit * std::conj(wave_dot(
                state.waves[q_index], state.velocity[p_index]));
        for (std::size_t component = 0; component < 3; ++component) {
            result[q_index][component] +=
                second_coefficient * target_cotangent[component];
        }
    }
    project_increment(result, state);
    return result;
}

SpectralIncrement SpectralDynamics::advection_jvp_fft(
    const SpectralState& state,
    const SpectralIncrement& direction) const {
    require_matching_increment(state, direction);
    return SpectralFftOperator::advection_jvp(
        state, direction, configuration_.compute_threads());
}

SpectralIncrement SpectralDynamics::advection_vjp_fft(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent) const {
    require_matching_increment(state, output_cotangent);
    SpectralIncrement cotangent = output_cotangent;
    project_increment(cotangent, state);
    SpectralIncrement result = SpectralFftOperator::advection_vjp(
        state, cotangent, configuration_.compute_threads());
    project_increment(result, state);
    return result;
}

SpectralIncrement SpectralDynamics::rhs(const SpectralState& state,
                                        SpectralReal viscosity) const {
    SpectralIncrement result = advection(state);
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[index]));
        for (std::size_t direction = 0; direction < 3; ++direction) {
            result[index][direction] =
                -result[index][direction] -
                viscosity * wave2 * state.velocity[index][direction];
        }
    }
    return result;
}

SpectralIncrement SpectralDynamics::rhs_jvp(
    const SpectralState& state, const SpectralIncrement& direction,
    SpectralReal viscosity) const {
    SpectralIncrement result = advection_jvp(state, direction);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                -result[mode][component] -
                viscosity * wave2 * direction[mode][component];
        }
    }
    return result;
}

SpectralIncrement SpectralDynamics::rhs_vjp(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent,
    SpectralReal viscosity) const {
    SpectralIncrement result = advection_vjp(state, output_cotangent);
    SpectralIncrement cotangent = output_cotangent;
    project_increment(cotangent, state);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                -result[mode][component] -
                viscosity * wave2 * cotangent[mode][component];
        }
    }
    project_increment(result, state);
    return result;
}

SpectralState SpectralDynamics::add_increment(
    const SpectralState& base, const SpectralIncrement& increment,
    SpectralReal factor) const {
    SpectralState result = base;
    for (std::size_t index = 0; index < result.waves.size(); ++index) {
        for (std::size_t direction = 0; direction < 3; ++direction) {
            result.velocity[index][direction] +=
                factor * increment[index][direction];
        }
    }
    return result;
}

void SpectralDynamics::enforce_constraints(SpectralState& state) const {
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        if (!is_positive_representative(wave)) {
            continue;
        }
        state.velocity[index] =
            project_divergence_free(wave, state.velocity[index]);
        state.velocity[state.index.at(-wave)] =
            conjugate(state.velocity[index]);
    }
}

void SpectralDynamics::project_increment(
    SpectralIncrement& increment, const SpectralState& layout) const {
    require_matching_increment(layout, increment);
    for (std::size_t mode = 0; mode < layout.waves.size(); ++mode) {
        const WaveVector wave = layout.waves[mode];
        if (!is_positive_representative(wave)) {
            continue;
        }
        const std::size_t negative = layout.index.at(-wave);
        ComplexVector combined{};
        for (std::size_t component = 0; component < 3; ++component) {
            combined[component] =
                0.5L * (increment[mode][component] +
                        std::conj(increment[negative][component]));
        }
        combined = project_divergence_free(wave, combined);
        increment[mode] = combined;
        increment[negative] = conjugate(combined);
    }
}

void SpectralDynamics::rk4_step(SpectralState& state, SpectralReal viscosity,
                                SpectralReal time_step) const {
    const SpectralIncrement k1 = rhs(state, viscosity);
    SpectralState stage = add_increment(state, k1, 0.5L * time_step);
    const SpectralIncrement k2 = rhs(stage, viscosity);
    stage = add_increment(state, k2, 0.5L * time_step);
    const SpectralIncrement k3 = rhs(stage, viscosity);
    stage = add_increment(state, k3, time_step);
    const SpectralIncrement k4 = rhs(stage, viscosity);

    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        for (std::size_t direction = 0; direction < 3; ++direction) {
            state.velocity[index][direction] +=
                (time_step / 6.0L) *
                (k1[index][direction] + 2.0L * k2[index][direction] +
                 2.0L * k3[index][direction] + k4[index][direction]);
        }
    }
    enforce_constraints(state);
}

SpectralIncrement SpectralDynamics::rk4_jvp(
    const SpectralState& state, const SpectralIncrement& direction,
    SpectralReal viscosity, SpectralReal time_step) const {
    require_matching_increment(state, direction);
    const SpectralIncrement k1 = rhs(state, viscosity);
    const SpectralIncrement dk1 = rhs_jvp(state, direction, viscosity);

    const SpectralState stage2 =
        add_increment(state, k1, 0.5L * time_step);
    SpectralIncrement direction2 = direction;
    add_scaled_increment(direction2, dk1, 0.5L * time_step);
    const SpectralIncrement k2 = rhs(stage2, viscosity);
    const SpectralIncrement dk2 = rhs_jvp(stage2, direction2, viscosity);

    const SpectralState stage3 =
        add_increment(state, k2, 0.5L * time_step);
    SpectralIncrement direction3 = direction;
    add_scaled_increment(direction3, dk2, 0.5L * time_step);
    const SpectralIncrement k3 = rhs(stage3, viscosity);
    const SpectralIncrement dk3 = rhs_jvp(stage3, direction3, viscosity);

    const SpectralState stage4 = add_increment(state, k3, time_step);
    SpectralIncrement direction4 = direction;
    add_scaled_increment(direction4, dk3, time_step);
    const SpectralIncrement dk4 = rhs_jvp(stage4, direction4, viscosity);

    SpectralIncrement result = direction;
    add_scaled_increment(result, dk1, time_step / 6.0L);
    add_scaled_increment(result, dk2, time_step / 3.0L);
    add_scaled_increment(result, dk3, time_step / 3.0L);
    add_scaled_increment(result, dk4, time_step / 6.0L);
    project_increment(result, state);
    return result;
}

SpectralIncrement SpectralDynamics::rk4_vjp(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent,
    SpectralReal viscosity, SpectralReal time_step) const {
    require_matching_increment(state, output_cotangent);
    const SpectralIncrement k1 = rhs(state, viscosity);
    const SpectralState stage2 =
        add_increment(state, k1, 0.5L * time_step);
    const SpectralIncrement k2 = rhs(stage2, viscosity);
    const SpectralState stage3 =
        add_increment(state, k2, 0.5L * time_step);
    const SpectralIncrement k3 = rhs(stage3, viscosity);
    const SpectralState stage4 = add_increment(state, k3, time_step);

    SpectralIncrement state_bar = output_cotangent;
    project_increment(state_bar, state);
    SpectralIncrement k1_bar =
        scaled_increment(state_bar, time_step / 6.0L);
    SpectralIncrement k2_bar =
        scaled_increment(state_bar, time_step / 3.0L);
    SpectralIncrement k3_bar =
        scaled_increment(state_bar, time_step / 3.0L);
    const SpectralIncrement k4_bar =
        scaled_increment(state_bar, time_step / 6.0L);

    const SpectralIncrement stage4_bar =
        rhs_vjp(stage4, k4_bar, viscosity);
    add_scaled_increment(state_bar, stage4_bar, 1.0L);
    add_scaled_increment(k3_bar, stage4_bar, time_step);

    const SpectralIncrement stage3_bar =
        rhs_vjp(stage3, k3_bar, viscosity);
    add_scaled_increment(state_bar, stage3_bar, 1.0L);
    add_scaled_increment(k2_bar, stage3_bar, 0.5L * time_step);

    const SpectralIncrement stage2_bar =
        rhs_vjp(stage2, k2_bar, viscosity);
    add_scaled_increment(state_bar, stage2_bar, 1.0L);
    add_scaled_increment(k1_bar, stage2_bar, 0.5L * time_step);

    const SpectralIncrement stage1_bar =
        rhs_vjp(state, k1_bar, viscosity);
    add_scaled_increment(state_bar, stage1_bar, 1.0L);
    project_increment(state_bar, state);
    return state_bar;
}

}  // namespace lemma
