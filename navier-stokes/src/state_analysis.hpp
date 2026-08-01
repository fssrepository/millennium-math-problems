#pragma once

#include "spectral_objective.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct StateAnalysisOptions {
    std::string state_path;
    std::string certificate_path;
    int top_modes = 12;
    SpectralReal active_relative_tolerance = 1e-24L;
};

struct StateShellAnalysis {
    int shell = 0;
    int modes = 0;
    int active_modes = 0;
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal q_gradient_norm2 = 0.0L;
};

struct StateModeAnalysis {
    WaveVector wave;
    SpectralReal energy = 0.0L;
    SpectralReal q_gradient_norm2 = 0.0L;
};

struct StateLineSample {
    SpectralReal step = 0.0L;
    SpectralReal q = 0.0L;
};

struct StateAnalysisReport {
    std::string state_path;
    int cutoff = 0;
    int modes = 0;
    int active_modes = 0;
    int highest_active_shell = 0;
    SpectralReal active_threshold = 0.0L;
    SpectralReal divergence_residual = 0.0L;
    SpectralReal reality_residual = 0.0L;
    SpectralReal projected_q_gradient_norm = 0.0L;
    SpectralReal retraction_directional_derivative = 0.0L;
    SpectralReal retraction_gradient_relative_error = 0.0L;
    SpectralReal homogeneous_h3_squared = 0.0L;
    SpectralReal homogeneous_h4_squared = 0.0L;
    StaticObjective objective;
    std::vector<StateShellAnalysis> shells;
    std::vector<StateModeAnalysis> top_modes;
    std::vector<StateLineSample> q_line_profile;
};

struct StateFamilyAnalysisOptions {
    std::string state_directory;
    std::string certificate_path;
    std::vector<int> cutoffs{1, 2, 3, 4};
};

struct StateFamilyAnalysisRow {
    int cutoff = 0;
    int modes = 0;
    SpectralReal energy = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal homogeneous_h3_squared = 0.0L;
    SpectralReal homogeneous_h4_squared = 0.0L;
    SpectralReal q = 0.0L;
    SpectralReal top_shell_energy = 0.0L;
    SpectralReal projection_residual = 0.0L;
};

struct StateFamilyAnalysisReport {
    std::string state_directory;
    SpectralReal top_shell_energy_exponent = 0.0L;
    SpectralReal maximum_projection_residual = 0.0L;
    std::vector<StateFamilyAnalysisRow> rows;
};

class SpectralStateReader {
public:
    [[nodiscard]] static SpectralState read_tsv(const std::string& path);
};

class SpectralStateAnalyzer {
public:
    [[nodiscard]] static StateAnalysisReport analyze(
        const SpectralState& state, const SpectralObjective& objective,
        const StateAnalysisOptions& options);
};

class StateAnalysisReporter {
public:
    static void write_console(const StateAnalysisReport& report,
                              std::ostream& out);
    static void write_json(const StateAnalysisReport& report,
                           std::ostream& out);
};

class StateFamilyAnalyzer {
public:
    [[nodiscard]] static StateFamilyAnalysisReport analyze(
        const StateFamilyAnalysisOptions& options,
        const SpectralObjective& objective);
};

class StateFamilyAnalysisReporter {
public:
    static void write_console(const StateFamilyAnalysisReport& report,
                              std::ostream& out);
    static void write_json(const StateFamilyAnalysisReport& report,
                           std::ostream& out);
};

int run_state_analysis(const StateAnalysisOptions& options, std::ostream& out);
int run_state_family_analysis(const StateFamilyAnalysisOptions& options,
                              std::ostream& out);

class StateAnalysisCli {
public:
    [[nodiscard]] static StateAnalysisOptions parse(int argc, char** argv,
                                                    int first);
    static void print_help(std::ostream& out);
};

class StateFamilyAnalysisCli {
public:
    [[nodiscard]] static StateFamilyAnalysisOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
};

}  // namespace lemma
