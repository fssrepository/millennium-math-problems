#include "local_quartic_shell_envelope.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace lemma {

LocalQuarticShellEnvelopeReport LocalQuarticShellEnvelope::analyze(
    const SpectralState& state,
    const LocalQuarticShellReport& shell_ledger) {
    LocalQuarticShellEnvelopeReport result;
    result.shells.resize(shell_ledger.dyadic_shells.size());
    for (std::size_t index = 0;
         index < shell_ledger.dyadic_shells.size(); ++index) {
        LocalQuarticShellEnvelopeRow& row = result.shells[index];
        row.shell = static_cast<int>(index);
        row.radius = std::ldexp(1.0L, row.shell);
        const std::size_t first = index == 0 ? 0 : index - 1;
        const std::size_t last = std::min(
            shell_ledger.dyadic_shells.size() - 1, index + 1);
        for (std::size_t neighbor = first; neighbor <= last; ++neighbor) {
            row.neighborhood_energy +=
                shell_ledger.dyadic_shells[neighbor].energy;
        }
        row.actual_local_advection_h1_squared = std::max(
            0.0L,
            shell_ledger.dyadic_shells[index]
                .local_advection_h1_squared);
        const SpectralReal radius2 = row.radius * row.radius;
        const SpectralReal radius4 = radius2 * radius2;
        const SpectralReal radius7 =
            radius4 * radius2 * row.radius;
        result.neighborhood_h3_moment +=
            radius2 * row.radius * row.neighborhood_energy;
        result.neighborhood_h4_moment +=
            radius4 * row.neighborhood_energy;
        result.shell_product_sum +=
            radius7 * row.neighborhood_energy *
            row.neighborhood_energy;
        result.actual_global_local_advection_h1_squared +=
            row.actual_local_advection_h1_squared;
        row.explicit_bound = result.explicit_constant * radius7 *
            row.neighborhood_energy * row.neighborhood_energy;
        if (row.explicit_bound > 0.0L) {
            row.bound_ratio =
                row.actual_local_advection_h1_squared /
                row.explicit_bound;
        }
        result.maximum_bound_ratio = std::max(
            result.maximum_bound_ratio, row.bound_ratio);
    }
    for (const LocalQuarticDyadicShellRow& row :
         shell_ledger.dyadic_shells) {
        result.state_enstrophy += row.enstrophy;
        result.state_palinstrophy += row.palinstrophy;
    }
    result.interpolated_h3_bound = std::sqrt(
        result.state_enstrophy * result.state_palinstrophy);
    result.global_shell_product_bound =
        result.explicit_constant * result.shell_product_sum;
    result.global_zp_bound = result.explicit_constant *
        result.h3_overlap_constant * result.h4_overlap_constant *
        result.interpolated_h3_bound * result.state_palinstrophy;
    if (result.global_shell_product_bound > 0.0L) {
        result.global_shell_product_ratio =
            result.actual_global_local_advection_h1_squared /
            result.global_shell_product_bound;
    }
    if (result.global_zp_bound > 0.0L) {
        result.global_zp_ratio =
            result.actual_global_local_advection_h1_squared /
            result.global_zp_bound;
    }
    std::vector<SpectralInteger> interaction_counts(state.waves.size(), 0);
    result.all_inputs_in_neighboring_shells = true;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (!TriadPartitioner::includes(
                state, interaction, TriadPartition::local)) {
            continue;
        }
        const auto [p_index, q_index, target_index] = interaction;
        ++interaction_counts[target_index];
        const SpectralInteger target_squared =
            norm_squared(state.waves[target_index]);
        const int shell = [&] {
            int value = 0;
            SpectralInteger next_squared = 4;
            while (target_squared >= next_squared) {
                ++value;
                next_squared *= 4;
            }
            return value;
        }();
        const SpectralReal radius = std::ldexp(1.0L, shell);
        const SpectralReal radius2 = radius * radius;
        result.maximum_target_frequency_ratio = std::max(
            result.maximum_target_frequency_ratio,
            static_cast<SpectralReal>(target_squared) /
                (4.0L * radius2));
        result.maximum_advected_frequency_ratio = std::max(
            result.maximum_advected_frequency_ratio,
            static_cast<SpectralReal>(
                norm_squared(state.waves[q_index])) /
                (16.0L * radius2));
        const auto input_is_neighbor = [&](std::size_t input_index) {
            const SpectralInteger input_squared =
                norm_squared(state.waves[input_index]);
            return input_squared * 4 >= radius2 &&
                   static_cast<SpectralReal>(input_squared) <
                       16.0L * radius2;
        };
        result.all_inputs_in_neighboring_shells =
            result.all_inputs_in_neighboring_shells &&
            input_is_neighbor(p_index) && input_is_neighbor(q_index);
    }
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralInteger target_squared =
            norm_squared(state.waves[mode]);
        int shell = 0;
        SpectralInteger next_squared = 4;
        while (target_squared >= next_squared) {
            ++shell;
            next_squared *= 4;
        }
        const SpectralReal radius = std::ldexp(1.0L, shell);
        const SpectralReal lattice_bound =
            std::pow(8.0L * radius + 1.0L, 3.0L);
        result.maximum_interaction_count_ratio = std::max(
            result.maximum_interaction_count_ratio,
            static_cast<SpectralReal>(interaction_counts[mode]) /
                lattice_bound);
    }
    result.all_geometry_checks_hold =
        result.all_inputs_in_neighboring_shells &&
        result.maximum_target_frequency_ratio <= 1.0L &&
        result.maximum_advected_frequency_ratio <= 1.0L &&
        result.maximum_interaction_count_ratio <= 1.0L;
    constexpr SpectralReal slack = 1.0e-15L;
    result.global_summation_holds =
        result.neighborhood_h3_moment <=
            result.h3_overlap_constant *
                result.interpolated_h3_bound * (1.0L + slack) &&
        result.neighborhood_h4_moment <=
            result.h4_overlap_constant * result.state_palinstrophy *
                (1.0L + slack) &&
        result.shell_product_sum <=
            result.neighborhood_h3_moment *
                result.neighborhood_h4_moment * (1.0L + slack) &&
        result.actual_global_local_advection_h1_squared <=
            result.global_shell_product_bound * (1.0L + slack) &&
        result.global_shell_product_bound <=
            result.global_zp_bound * (1.0L + slack);
    result.all_bounds_hold = result.cutoff_independent &&
        result.all_geometry_checks_hold &&
        result.global_summation_holds &&
        result.maximum_bound_ratio <= 1.0L + slack;
    return result;
}

}  // namespace lemma
