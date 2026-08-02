#include "adversary_reporter.hpp"

#include <iomanip>
#include <ostream>

namespace lemma {

void AdversaryReporter::write_console(const AdversaryReport& report,
                                      std::ostream& out) {
    out << "# workers=" << report.workers << ", backend=" << report.backend << '\n'
        << "cutoff,modes,evaluations,accepted,E,Z,P,V,D,Q=D^4Z,I=D^4Z^2\n"
        << std::setprecision(12);
    for (const auto& row : report.rows) {
        out << row.cutoff << ',' << row.modes << ',' << row.evaluations << ','
            << row.accepted_mutations << ',' << static_cast<double>(row.energy) << ','
            << static_cast<double>(row.enstrophy) << ','
            << static_cast<double>(row.palinstrophy) << ','
            << static_cast<double>(row.vortex_stretching) << ','
            << static_cast<double>(row.depletion) << ','
            << static_cast<double>(row.q) << ','
            << static_cast<double>(row.critical_integrand) << '\n';
    }
    out << "\nQ cutoff growth ratio: " << static_cast<double>(report.q_growth_ratio)
        << "\nQ fitted cutoff exponent: "
        << static_cast<double>(report.q_cutoff_log_slope)
        << "\nWarm-start embedding monotonicity: "
        << (report.embedding_monotonicity ? "PASS" : "FAIL")
        << "\nInterpretation: positive finite-cutoff growth attacks the strong "
           "pointwise quarter-depletion lemma. It does not by itself falsify "
           "the time-integrated L4-A statement.\n"
        << "\nDynamic Galerkin check (nu=" << static_cast<double>(report.viscosity)
        << ", T=" << static_cast<double>(report.time)
        << ", objective=" << report.dynamic_objective
        << ", optimizer=" << report.dynamic_optimizer
        << ", gradient_method=" << report.gradient_method
        << ", minimum_dyadic_gap=" << report.minimum_dyadic_gap
        << ", H" << report.sobolev_order << "_cap="
        << static_cast<double>(report.sobolev_cap)
        << ", critical_density_shift="
        << static_cast<double>(report.critical_density_shift)
        << ", dynamic_restarts=" << report.dynamic_restarts
        << ", replay_each_cutoff="
        << (report.dynamic_replay_each_cutoff ? "yes" : "no") << ")\n"
        << "cutoff,steps,int_D4Z2_refined,int_local_D4Z2,int_nonlocal_D4Z2,"
           "int_near_nonlocal_D4Z2,int_far_nonlocal_D4Z2,"
           "int_selected_gap_tail_D4Z2,"
           "dt_relative_error,search_obj_initial,search_obj_final,"
           "initial_local_density,final_local_density,initial_Z,initial_P,E0P0,"
           "local_log_gain,local_log_gain_over_Tk0Z,"
           "EP_shifted_log_gain,EP_shifted_log_gain_over_Tk0Z,"
           "initial_D4Z,final_D4Z,log_Q_gain,max_D4Z,max_local_D4Z,max_nonlocal_D4Z,"
           "max_positive_dlogQ_over_k0Z,q_derivative_error,strong_L4_envelope,"
           "envelope_use,max_Z,max_omega_inf,max_holder_half,"
           "max_stretch_alignment,nonlocal_V_fraction,V_partition_residual,"
           "final_E,energy_balance_residual\n";
    for (const auto& row : report.rows) {
        out << row.cutoff << ',' << row.dynamic_steps << ','
            << static_cast<double>(row.dynamic_integral) << ','
            << static_cast<double>(row.dynamic_local_integral) << ','
            << static_cast<double>(row.dynamic_nonlocal_integral) << ','
            << static_cast<double>(row.dynamic_near_nonlocal_integral) << ','
            << static_cast<double>(row.dynamic_far_nonlocal_integral) << ','
            << static_cast<double>(row.dynamic_selected_gap_tail_integral)
            << ','
            << static_cast<double>(row.dynamic_dt_relative_error) << ','
            << static_cast<double>(row.dynamic_search_initial_objective) << ','
            << static_cast<double>(row.dynamic_search_final_objective) << ','
            << static_cast<double>(
                   row.dynamic_initial_local_critical_density) << ','
            << static_cast<double>(
                   row.dynamic_final_local_critical_density) << ','
            << static_cast<double>(row.dynamic_initial_enstrophy) << ','
            << static_cast<double>(row.dynamic_initial_palinstrophy) << ','
            << static_cast<double>(row.dynamic_initial_ep_shift) << ','
            << static_cast<double>(
                   row.dynamic_local_critical_log_gain) << ','
            << static_cast<double>(
                   row.dynamic_local_log_gain_rate_ratio) << ','
            << static_cast<double>(
                   row.dynamic_ep_shifted_local_log_gain) << ','
            << static_cast<double>(
                   row.dynamic_ep_shifted_log_gain_rate_ratio) << ','
            << static_cast<double>(row.dynamic_initial_q) << ','
            << static_cast<double>(row.dynamic_final_q) << ','
            << static_cast<double>(row.dynamic_log_q_gain) << ','
            << static_cast<double>(row.dynamic_maximum_q) << ','
            << static_cast<double>(row.dynamic_maximum_local_q) << ','
            << static_cast<double>(row.dynamic_maximum_nonlocal_q) << ','
            << static_cast<double>(row.dynamic_q_log_growth_ratio) << ','
            << static_cast<double>(row.dynamic_q_derivative_error) << ','
            << static_cast<double>(row.strong_l4_envelope) << ','
            << static_cast<double>(row.envelope_utilization) << ','
            << static_cast<double>(row.dynamic_maximum_enstrophy) << ','
            << static_cast<double>(row.dynamic_maximum_vorticity) << ','
            << static_cast<double>(row.dynamic_maximum_holder_half) << ','
            << static_cast<double>(row.dynamic_maximum_stretch_alignment) << ','
            << static_cast<double>(row.dynamic_nonlocal_vortex_fraction) << ','
            << static_cast<double>(row.dynamic_partition_residual) << ','
            << static_cast<double>(row.dynamic_final_energy) << ','
            << static_cast<double>(row.dynamic_energy_balance_residual)
            << "  # dynamic_evals=" << row.dynamic_evaluations
            << ", winning_restart=" << row.dynamic_winning_restart
            << ", accepted_mutations=" << row.dynamic_accepted_mutations
            << ", accepted_gradient="
            << row.dynamic_accepted_gradient_steps
            << ", final_sobolev="
            << static_cast<double>(row.dynamic_sobolev_value);
        if (!row.dynamic_gradient_trace.empty()) {
            const auto& final_trace = row.dynamic_gradient_trace.back();
            out << ", final_projected_gradient="
                << static_cast<double>(
                       final_trace.projected_gradient_norm)
                << ", final_step="
                << static_cast<double>(final_trace.accepted_step);
        }
        out << '\n';
    }
}

void AdversaryReporter::write_json(const AdversaryReport& report,
                                   std::ostream& out) {
    out << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-l4-adversary-v2\",\n"
        << "  \"candidate\": \"cutoff-uniform trajectory bound on Q=D^4 Z\",\n"
        << "  \"optimizer\": {\"restarts\": " << report.restarts
        << ", \"dynamic_restarts\": " << report.dynamic_restarts
        << ", \"dynamic_replay_each_cutoff\": "
        << (report.dynamic_replay_each_cutoff ? "true" : "false")
        << ", \"generations\": " << report.generations
        << ", \"dynamic_generations\": " << report.dynamic_generations
        << ", \"dynamic_objective\": \"" << report.dynamic_objective << "\""
        << ", \"dynamic_optimizer\": \"" << report.dynamic_optimizer << "\""
        << ", \"gradient_method\": \"" << report.gradient_method << "\""
        << ", \"minimum_dyadic_gap\": " << report.minimum_dyadic_gap
        << ", \"sobolev_order\": " << report.sobolev_order
        << ", \"sobolev_cap\": "
        << static_cast<double>(report.sobolev_cap)
        << ", \"critical_density_shift\": "
        << static_cast<double>(report.critical_density_shift)
        << ", \"mutation\": " << static_cast<double>(report.mutation)
        << ", \"seed\": " << report.seed << "},\n"
        << "  \"threads\": " << report.workers << ",\n"
        << "  \"backend\": \"" << report.backend << "\",\n"
        << "  \"evolution\": {\"viscosity\": "
        << static_cast<double>(report.viscosity)
        << ", \"time\": " << static_cast<double>(report.time)
        << ", \"requested_dt\": " << static_cast<double>(report.requested_dt)
        << "},\n  \"results\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const auto& row = report.rows[index];
        out << "    {\"cutoff\": " << row.cutoff
            << ", \"modes\": " << row.modes
            << ", \"evaluations\": " << row.evaluations
            << ", \"accepted_mutations\": " << row.accepted_mutations
            << ", \"E\": " << static_cast<double>(row.energy)
            << ", \"Z\": " << static_cast<double>(row.enstrophy)
            << ", \"P\": " << static_cast<double>(row.palinstrophy)
            << ", \"V\": " << static_cast<double>(row.vortex_stretching)
            << ", \"D\": " << static_cast<double>(row.depletion)
            << ", \"Q\": " << static_cast<double>(row.q)
            << ", \"critical_integrand\": "
            << static_cast<double>(row.critical_integrand)
            << ", \"dynamic_integral_D4Z2\": "
            << static_cast<double>(row.dynamic_integral)
            << ", \"dynamic_coarse_integral_D4Z2\": "
            << static_cast<double>(row.dynamic_coarse_integral)
            << ", \"dynamic_dt_relative_error\": "
            << static_cast<double>(row.dynamic_dt_relative_error)
            << ", \"dynamic_search_initial_objective\": "
            << static_cast<double>(row.dynamic_search_initial_objective)
            << ", \"dynamic_search_final_objective\": "
            << static_cast<double>(row.dynamic_search_final_objective)
            << ", \"dynamic_initial_local_critical_density\": "
            << static_cast<double>(
                   row.dynamic_initial_local_critical_density)
            << ", \"dynamic_final_local_critical_density\": "
            << static_cast<double>(
                   row.dynamic_final_local_critical_density)
            << ", \"dynamic_initial_Z\": "
            << static_cast<double>(row.dynamic_initial_enstrophy)
            << ", \"dynamic_initial_P\": "
            << static_cast<double>(row.dynamic_initial_palinstrophy)
            << ", \"dynamic_initial_E0P0_shift\": "
            << static_cast<double>(row.dynamic_initial_ep_shift)
            << ", \"dynamic_local_critical_log_gain\": "
            << static_cast<double>(row.dynamic_local_critical_log_gain)
            << ", \"dynamic_local_log_gain_over_Tk0Z\": "
            << static_cast<double>(row.dynamic_local_log_gain_rate_ratio)
            << ", \"dynamic_E0P0_shifted_local_log_gain\": "
            << static_cast<double>(row.dynamic_ep_shifted_local_log_gain)
            << ", \"dynamic_E0P0_shifted_log_gain_over_Tk0Z\": "
            << static_cast<double>(
                   row.dynamic_ep_shifted_log_gain_rate_ratio)
            << ", \"dynamic_max_D4Z\": "
            << static_cast<double>(row.dynamic_maximum_q)
            << ", \"dynamic_initial_D4Z\": "
            << static_cast<double>(row.dynamic_initial_q)
            << ", \"dynamic_final_D4Z\": "
            << static_cast<double>(row.dynamic_final_q)
            << ", \"dynamic_log_Q_gain\": "
            << static_cast<double>(row.dynamic_log_q_gain)
            << ", \"dynamic_max_local_D4Z\": "
            << static_cast<double>(row.dynamic_maximum_local_q)
            << ", \"dynamic_max_nonlocal_D4Z\": "
            << static_cast<double>(row.dynamic_maximum_nonlocal_q)
            << ", \"dynamic_max_near_nonlocal_D4Z\": "
            << static_cast<double>(row.dynamic_maximum_near_nonlocal_q)
            << ", \"dynamic_max_far_nonlocal_D4Z\": "
            << static_cast<double>(row.dynamic_maximum_far_nonlocal_q)
            << ", \"dynamic_max_positive_dlogQ_over_k0Z\": "
            << static_cast<double>(row.dynamic_q_log_growth_ratio)
            << ", \"dynamic_q_derivative_refinement_error\": "
            << static_cast<double>(row.dynamic_q_derivative_error)
            << ", \"dynamic_max_Z\": "
            << static_cast<double>(row.dynamic_maximum_enstrophy)
            << ", \"dynamic_max_omega_inf\": "
            << static_cast<double>(row.dynamic_maximum_vorticity)
            << ", \"dynamic_max_holder_half_coherence\": "
            << static_cast<double>(row.dynamic_maximum_holder_half)
            << ", \"dynamic_max_stretch_alignment\": "
            << static_cast<double>(row.dynamic_maximum_stretch_alignment)
            << ", \"dynamic_geometry_samples\": " << row.dynamic_geometry_samples
            << ", \"dynamic_integral_abs_local_V\": "
            << static_cast<double>(row.dynamic_integral_absolute_local_vortex)
            << ", \"dynamic_integral_abs_nonlocal_V\": "
            << static_cast<double>(row.dynamic_integral_absolute_nonlocal_vortex)
            << ", \"dynamic_integral_abs_total_V\": "
            << static_cast<double>(row.dynamic_integral_absolute_total_vortex)
            << ", \"dynamic_integral_local_D4Z2\": "
            << static_cast<double>(row.dynamic_local_integral)
            << ", \"dynamic_integral_nonlocal_D4Z2\": "
            << static_cast<double>(row.dynamic_nonlocal_integral)
            << ", \"dynamic_integral_near_nonlocal_D4Z2\": "
            << static_cast<double>(row.dynamic_near_nonlocal_integral)
            << ", \"dynamic_integral_far_nonlocal_D4Z2\": "
            << static_cast<double>(row.dynamic_far_nonlocal_integral)
            << ", \"dynamic_integral_selected_gap_tail_D4Z2\": "
            << static_cast<double>(row.dynamic_selected_gap_tail_integral)
            << ", \"dynamic_max_selected_gap_tail_D4Z\": "
            << static_cast<double>(row.dynamic_maximum_selected_gap_tail_q)
            << ", \"dynamic_max_V_partition_residual\": "
            << static_cast<double>(row.dynamic_partition_residual)
            << ", \"dynamic_final_E\": "
            << static_cast<double>(row.dynamic_final_energy)
            << ", \"dynamic_energy_balance_residual\": "
            << static_cast<double>(row.dynamic_energy_balance_residual)
            << ", \"dynamic_evaluations\": " << row.dynamic_evaluations
            << ", \"dynamic_winning_restart\": "
            << row.dynamic_winning_restart
            << ", \"dynamic_accepted_mutations\": "
            << row.dynamic_accepted_mutations
            << ", \"dynamic_accepted_gradient_steps\": "
            << row.dynamic_accepted_gradient_steps
            << ", \"dynamic_final_sobolev_value\": "
            << static_cast<double>(row.dynamic_sobolev_value)
            << ", \"dynamic_restart_objectives\": [";
        for (std::size_t restart = 0;
             restart < row.dynamic_restart_objectives.size(); ++restart) {
            out << static_cast<double>(
                       row.dynamic_restart_objectives[restart])
                << (restart + 1 == row.dynamic_restart_objectives.size()
                        ? ""
                        : ", ");
        }
        out << "]"
            << ", \"dynamic_gradient_trace\": [";
        for (std::size_t trace_index = 0;
             trace_index < row.dynamic_gradient_trace.size();
             ++trace_index) {
            const auto& point = row.dynamic_gradient_trace[trace_index];
            out << "{\"iteration\": " << point.iteration
                << ", \"objective_before\": "
                << static_cast<double>(point.objective_before)
                << ", \"objective_after\": "
                << static_cast<double>(point.objective_after)
                << ", \"projected_gradient_norm\": "
                << static_cast<double>(point.projected_gradient_norm)
                << ", \"accepted_step\": "
                << static_cast<double>(point.accepted_step)
                << ", \"sobolev_value\": "
                << static_cast<double>(point.sobolev_value)
                << ", \"line_search_evaluations\": "
                << point.line_search_evaluations
                << ", \"used_steepest_fallback\": "
                << (point.used_steepest_fallback ? "true" : "false")
                << ", \"accepted\": "
                << (point.accepted ? "true" : "false") << '}'
                << (trace_index + 1 == row.dynamic_gradient_trace.size()
                        ? ""
                        : ", ");
        }
        out << "]}"
            << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    out << "  ],\n"
        << "  \"Q_growth_ratio\": " << static_cast<double>(report.q_growth_ratio)
        << ",\n  \"Q_cutoff_log_slope\": "
        << static_cast<double>(report.q_cutoff_log_slope)
        << ",\n  \"embedding_monotonicity\": "
        << (report.embedding_monotonicity ? "true" : "false") << ",\n"
        << "  \"logical_status\": \"finite adversarial evidence only; derive an analytic family before rejection\"\n"
        << "}\n";
}

}  // namespace lemma
