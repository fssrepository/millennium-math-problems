#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSignatureTrajectoryOptions {
    std::string state_path;
    std::string certificate_path;
    int minimum_cutoff = 3;
    int maximum_cutoff = 6;
    int steps = 10;
    int time_refinement = 2;
    int workers = 12;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.001L;
};

struct LocalSignatureTrajectoryRow {
    int cutoff = 0;
    SpectralReal initial_amplification = 0.0L;
    SpectralReal maximum_amplification = 0.0L;
    SpectralReal final_amplification = 0.0L;
    SpectralReal critical_integral = 0.0L;
    SpectralReal refined_critical_integral = 0.0L;
    SpectralReal time_refinement_relative_difference = 0.0L;
    SpectralReal square_signature_integral = 0.0L;
    SpectralReal maximum_critical_density = 0.0L;
    SpectralReal maximum_square_signature_density = 0.0L;
    SpectralReal log_factor_correlation = 0.0L;
    SpectralReal square_density_at_amplification_peak = 0.0L;
    SpectralReal critical_density_at_amplification_peak = 0.0L;
    SpectralReal amplification_at_square_density_peak = 0.0L;
    SpectralReal critical_density_at_square_density_peak = 0.0L;
    SpectralReal amplification_at_critical_peak = 0.0L;
    SpectralReal square_density_at_critical_peak = 0.0L;
    int amplification_peak_step = 0;
    int square_density_peak_step = 0;
    int critical_peak_step = 0;
    SpectralReal maximum_factorization_residual = 0.0L;
};

struct LocalSignatureTrajectoryReport {
    std::string state_path;
    int workers = 0;
    int kernel_threads = 1;
    int steps = 0;
    int time_refinement = 1;
    SpectralReal viscosity = 0.0L;
    SpectralReal time_step = 0.0L;
    SpectralReal last_relative_critical_integral_difference = 0.0L;
    std::vector<LocalSignatureTrajectoryRow> rows;
};

class LocalSignatureTrajectoryAnalyzer {
public:
    [[nodiscard]] static LocalSignatureTrajectoryReport analyze(
        const LocalSignatureTrajectoryOptions& options);
};

class LocalSignatureTrajectoryCli {
public:
    [[nodiscard]] static LocalSignatureTrajectoryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const LocalSignatureTrajectoryOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
