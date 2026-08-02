#include "family_reporter.hpp"

#include <iomanip>
#include <ostream>

namespace lemma {

void FamilyReporter::write_console(const FamilyReport& report, std::ostream& out) {
    out << "# workers=" << report.threads << ", backend=" << report.backend << '\n'
        << "seed,cutoff,modes,E0,Z0,int_D4Z2,max_D4Z,max_local_D4Z,"
           "max_nonlocal_D4Z,max_positive_dlogQ_over_k0Z,q_derivative_error,"
           "strong_L4_envelope,envelope_use,dt_abs_error,dt_rel_error,max_Z,"
           "max_omega_inf,max_holder_half,local_integral,nonlocal_integral,"
           "projection_residual\n";
    out << std::setprecision(12);
    for (const auto& row : report.runs) {
        out << row.seed << ',' << row.cutoff << ',' << row.modes << ','
            << static_cast<double>(row.initial_energy) << ','
            << static_cast<double>(row.initial_enstrophy) << ','
            << static_cast<double>(row.integral_critical) << ','
            << static_cast<double>(row.maximum_q) << ','
            << static_cast<double>(row.maximum_local_q) << ','
            << static_cast<double>(row.maximum_nonlocal_q) << ','
            << static_cast<double>(row.maximum_positive_q_log_growth_ratio) << ','
            << static_cast<double>(row.q_derivative_refinement_error) << ','
            << static_cast<double>(row.energy_identity_envelope) << ','
            << static_cast<double>(row.envelope_utilization) << ','
            << static_cast<double>(row.dt_absolute_error) << ','
            << static_cast<double>(row.dt_relative_error) << ','
            << static_cast<double>(row.maximum_enstrophy) << ','
            << static_cast<double>(row.maximum_vorticity) << ','
            << static_cast<double>(row.maximum_holder_half) << ','
            << static_cast<double>(row.local_integral) << ','
            << static_cast<double>(row.nonlocal_integral) << ','
            << static_cast<double>(row.projection_residual) << '\n';
    }
    out << "\nseed,last_integral_relative_increment,max_Q_over_cutoffs,"
           "endpoint_Q_growth,endpoint_Q_exponent,tail_record_growth,"
           "tail_record_exponent\n";
    for (const auto& row : report.summaries) {
        out << row.seed << ','
            << static_cast<double>(row.last_relative_increment) << ','
            << static_cast<double>(row.maximum_q) << ','
            << static_cast<double>(row.endpoint_q_growth_ratio) << ','
            << static_cast<double>(row.endpoint_q_log_slope) << ','
            << static_cast<double>(row.tail_record_growth_ratio) << ','
            << static_cast<double>(row.tail_record_log_slope) << '\n';
    }
    out << "worst max-Q tail-record exponent: "
        << static_cast<double>(report.worst_tail_log_slope)
        << " (seed " << report.worst_tail_seed << ")\n";
}

void FamilyReporter::write_json(const FamilyReport& report, std::ostream& out) {
    out << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-projective-family-v2\",\n"
        << "  \"initial_family\": \"deterministic divergence-free analytic Fourier field\",\n"
        << "  \"seed\": " << report.initial_seed << ",\n"
        << "  \"seed_count\": " << report.seed_count << ",\n"
        << "  \"spectral_decay\": " << static_cast<double>(report.spectral_decay) << ",\n"
        << "  \"viscosity\": " << static_cast<double>(report.viscosity) << ",\n"
        << "  \"time\": " << static_cast<double>(report.time) << ",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"backend\": \"" << report.backend << "\",\n"
        << "  \"runs\": [\n";
    for (std::size_t index = 0; index < report.runs.size(); ++index) {
        const auto& row = report.runs[index];
        out << "    {\"seed\": " << row.seed
            << ", \"cutoff\": " << row.cutoff
            << ", \"modes\": " << row.modes
            << ", \"E0\": " << static_cast<double>(row.initial_energy)
            << ", \"Z0\": " << static_cast<double>(row.initial_enstrophy)
            << ", \"integral_D4Z2\": " << static_cast<double>(row.integral_critical)
            << ", \"max_D4Z\": " << static_cast<double>(row.maximum_q)
            << ", \"max_local_D4Z\": " << static_cast<double>(row.maximum_local_q)
            << ", \"max_nonlocal_D4Z\": " << static_cast<double>(row.maximum_nonlocal_q)
            << ", \"max_positive_dlogQ_over_k0Z\": "
            << static_cast<double>(row.maximum_positive_q_log_growth_ratio)
            << ", \"q_derivative_refinement_error\": "
            << static_cast<double>(row.q_derivative_refinement_error)
            << ", \"q_derivative_samples\": " << row.q_derivative_samples
            << ", \"q_enstrophy_envelope\": "
            << static_cast<double>(row.q_enstrophy_envelope)
            << ", \"energy_identity_envelope\": "
            << static_cast<double>(row.energy_identity_envelope)
            << ", \"envelope_utilization\": "
            << static_cast<double>(row.envelope_utilization)
            << ", \"factorization_violation\": "
            << static_cast<double>(row.factorization_violation)
            << ", \"dt_relative_error\": "
            << static_cast<double>(row.dt_relative_error)
            << ", \"dt_absolute_error\": "
            << static_cast<double>(row.dt_absolute_error)
            << ", \"max_Z\": " << static_cast<double>(row.maximum_enstrophy)
            << ", \"max_omega_inf\": " << static_cast<double>(row.maximum_vorticity)
            << ", \"max_holder_half\": " << static_cast<double>(row.maximum_holder_half)
            << ", \"local_integral\": " << static_cast<double>(row.local_integral)
            << ", \"nonlocal_integral\": " << static_cast<double>(row.nonlocal_integral)
            << ", \"projection_residual\": "
            << static_cast<double>(row.projection_residual) << '}'
            << (index + 1 == report.runs.size() ? "\n" : ",\n");
    }
    out << "  ],\n  \"family_summaries\": [\n";
    for (std::size_t index = 0; index < report.summaries.size(); ++index) {
        const auto& row = report.summaries[index];
        out << "    {\"seed\": " << row.seed
            << ", \"last_cutoff_increment\": "
            << static_cast<double>(row.last_increment)
            << ", \"last_cutoff_relative_increment\": "
            << static_cast<double>(row.last_relative_increment)
            << ", \"max_Q_growth_ratio\": "
            << static_cast<double>(row.endpoint_q_growth_ratio)
            << ", \"max_Q_cutoff_log_slope\": "
            << static_cast<double>(row.endpoint_q_log_slope)
            << ", \"max_Q_over_cutoffs\": " << static_cast<double>(row.maximum_q)
            << ", \"tail_record_growth_ratio\": "
            << static_cast<double>(row.tail_record_growth_ratio)
            << ", \"tail_record_log_slope\": "
            << static_cast<double>(row.tail_record_log_slope) << '}'
            << (index + 1 == report.summaries.size() ? "\n" : ",\n");
    }
    out << "  ],\n"
        << "  \"worst_last_cutoff_relative_increment\": "
        << static_cast<double>(report.worst_last_relative_increment) << ",\n"
        << "  \"worst_tail_record_growth_ratio\": "
        << static_cast<double>(report.worst_tail_growth_ratio) << ",\n"
        << "  \"worst_tail_record_log_slope\": "
        << static_cast<double>(report.worst_tail_log_slope) << ",\n"
        << "  \"strong_L4_reduction\": \"sup_N,t Q_N times E0/(2 nu) bounds integral D_N^4 Z_N^2 dt\",\n"
        << "  \"logical_status\": \"finite projective convergence evidence; not an infinite-cutoff proof\"\n"
        << "}\n";
}

}  // namespace lemma
