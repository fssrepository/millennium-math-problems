#pragma once

#include "spectral_state.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSldProjectiveNormalizationAlternatingPhase {
    int cycle = 0;
    std::string component;
    SpectralReal component_initial = 0.0L;
    SpectralReal component_final = 0.0L;
    SpectralReal open_initial = 0.0L;
    SpectralReal open_final = 0.0L;
    SpectralReal component_projected_gradient_norm = 0.0L;
    SpectralReal open_projected_gradient_norm = 0.0L;
    int component_accepted_steps = 0;
    int open_accepted_steps = 0;
    int component_evaluations = 0;
    int open_evaluations = 0;
    std::string component_state_path;
    std::string open_state_path;
};

struct LocalSldProjectiveNormalizationAlternatingReport {
    int cutoff = 0;
    SpectralInteger core_maximum_height = 0;
    std::string component;
    SpectralReal initial_open_power_one = 0.0L;
    SpectralReal final_open_power_one = 0.0L;
    SpectralReal improvement_factor = 0.0L;
    SpectralState state;
    std::vector<LocalSldProjectiveNormalizationAlternatingPhase> phases;
    bool every_phase_finite = false;
    bool finite_alternating_search_is_not_a_proof = true;
    bool candidate_lemma_proved = false;
};

struct LocalSldProjectiveNormalizationAlternatingOptions {
    std::string state_path;
    std::string output_state_path;
    std::string certificate_path;
    std::string selection = "double-triple-remainder-without-123";
    std::string component = "dominant";
    SpectralInteger core_maximum_height = 8;
    int cycles = 2;
    int component_iterations = 12;
    int open_iterations = 16;
    int line_search_steps = 12;
    int lbfgs_history = 8;
    int threads = 12;
    SpectralReal initial_step = 0.1L;
};

class LocalSldProjectiveNormalizationAlternatingAdversaryCli {
public:
    [[nodiscard]] static
    LocalSldProjectiveNormalizationAlternatingOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveNormalizationAlternatingOptions& options,
        std::ostream& out);
};

}  // namespace lemma
