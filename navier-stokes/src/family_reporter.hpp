#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct FamilyReportRow {
    std::uint64_t seed = 0;
    int cutoff = 0;
    std::size_t modes = 0;
    long double initial_energy = 0.0L;
    long double initial_enstrophy = 0.0L;
    long double integral_critical = 0.0L;
    long double maximum_q = 0.0L;
    long double maximum_local_q = 0.0L;
    long double maximum_nonlocal_q = 0.0L;
    long double maximum_positive_q_log_growth_ratio = 0.0L;
    long double q_derivative_refinement_error = 0.0L;
    int q_derivative_samples = 0;
    long double q_enstrophy_envelope = 0.0L;
    long double energy_identity_envelope = 0.0L;
    long double envelope_utilization = 0.0L;
    long double factorization_violation = 0.0L;
    long double dt_absolute_error = 0.0L;
    long double dt_relative_error = 0.0L;
    long double maximum_enstrophy = 0.0L;
    long double maximum_vorticity = 0.0L;
    long double maximum_holder_half = 0.0L;
    long double local_integral = 0.0L;
    long double nonlocal_integral = 0.0L;
    long double projection_residual = 0.0L;
};

struct FamilySummaryRow {
    std::uint64_t seed = 0;
    long double last_increment = 0.0L;
    long double last_relative_increment = 0.0L;
    long double endpoint_q_growth_ratio = 1.0L;
    long double endpoint_q_log_slope = 0.0L;
    long double maximum_q = 0.0L;
    long double tail_record_growth_ratio = 1.0L;
    long double tail_record_log_slope = 0.0L;
};

struct FamilyReport {
    std::uint64_t initial_seed = 0;
    int seed_count = 1;
    long double spectral_decay = 0.0L;
    long double viscosity = 0.0L;
    long double time = 0.0L;
    int threads = 1;
    std::string backend;
    std::vector<FamilyReportRow> runs;
    std::vector<FamilySummaryRow> summaries;
    std::uint64_t worst_tail_seed = 0;
    long double worst_tail_growth_ratio = 1.0L;
    long double worst_tail_log_slope = 0.0L;
    long double worst_last_relative_increment = 0.0L;
};

class FamilyReporter {
public:
    static void write_console(const FamilyReport& report, std::ostream& out);
    static void write_json(const FamilyReport& report, std::ostream& out);
};

}  // namespace lemma
