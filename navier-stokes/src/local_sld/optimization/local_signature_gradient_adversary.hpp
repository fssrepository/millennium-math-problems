#pragma once

#include "local_signature_state_factory.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSignatureGradientOptions {
    std::string objective = "amplification";
    int iterations = 8;
    int line_search_steps = 14;
    SpectralReal initial_step = 0.1L;
};

struct LocalSignatureGradientTraceRow {
    int iteration = 0;
    SpectralReal objective = 0.0L;
    SpectralReal projected_gradient_norm = 0.0L;
    SpectralReal accepted_step = 0.0L;
    bool accepted = false;
};

struct LocalSignatureGradientResult {
    SpectralState state;
    SpectralReal initial_objective = 0.0L;
    SpectralReal objective = 0.0L;
    int accepted_steps = 0;
    int evaluations = 0;
    std::vector<LocalSignatureGradientTraceRow> trace;
};

class LocalSignatureGradientAdversary {
public:
    [[nodiscard]] static LocalSignatureGradientResult maximize(
        const SpectralState& initial,
        const LocalSignatureGradientOptions& options);
};

struct LocalSignatureGradientCliOptions {
    int minimum_cutoff = 3;
    int maximum_cutoff = 6;
    int restarts = 12;
    int workers = 12;
    int iterations = 8;
    int line_search_steps = 14;
    SpectralReal initial_step = 0.1L;
    std::uint64_t seed = 20260801;
    std::string profile = "flat";
    std::string objective = "amplification";
    std::string certificate_path;
};

struct LocalSignatureGradientCutoffRow {
    int cutoff = 0;
    SpectralReal best_initial_objective = 0.0L;
    SpectralReal best_objective = 0.0L;
    SpectralReal improvement_factor = 0.0L;
    SpectralReal absolute_signed_transfer = 0.0L;
    SpectralReal signature_transfer_l2 = 0.0L;
    SpectralReal normalized_lsf2_l2 = 0.0L;
    SpectralReal normalized_viscous_transfer = 0.0L;
    std::uint64_t best_seed = 0;
    int accepted_steps = 0;
    int evaluations = 0;
};

struct LocalSignatureGradientReport {
    std::string profile;
    std::string objective;
    int workers = 0;
    int restarts = 0;
    int iterations = 0;
    SpectralReal fitted_log_slope = 0.0L;
    SpectralReal critical_log_slope = 0.5L;
    bool finite_search_is_not_a_proof = true;
    std::vector<LocalSignatureGradientCutoffRow> rows;
};

class LocalSignatureGradientCli {
public:
    [[nodiscard]] static LocalSignatureGradientCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const LocalSignatureGradientCliOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
