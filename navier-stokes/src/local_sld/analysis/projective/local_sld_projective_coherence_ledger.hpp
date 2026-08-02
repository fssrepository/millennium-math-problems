#pragma once

#include "spectral_dynamics.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>

namespace lemma {

struct LocalSldProjectiveCoherenceReport {
    int cutoff = 0;
    int threads = 1;
    std::size_t selected_interactions = 0;
    std::size_t projective_shape_count = 0;
    std::size_t maximum_output_shape_multiplicity = 0;
    bool excludes_signature_123 = false;
    bool excludes_triple_family = false;
    SpectralReal total_advection_norm2 = 0.0L;
    SpectralReal projective_square_function_norm2 = 0.0L;
    SpectralReal coherent_synthesis_ratio = 0.0L;
    SpectralReal coherent_synthesis_amplification = 0.0L;
    SpectralReal maximum_output_synthesis_ratio = 0.0L;
    SpectralReal maximum_output_synthesis_amplification = 0.0L;
    SpectralReal maximum_output_coherent_fraction = 0.0L;
    SpectralReal maximum_output_square_function_fraction = 0.0L;
    WaveVector maximum_output_wave{};
    SpectralReal effective_advection_shapes = 0.0L;
    SpectralReal dominant_advection_shape_fraction = 0.0L;
    SpectralReal reconstruction_relative_error = 0.0L;
    bool exact_projective_reconstruction = false;
    bool cutoff_independent_synthesis_bound_proved = false;
};

class LocalSldProjectiveCoherenceLedger {
public:
    [[nodiscard]] static LocalSldProjectiveCoherenceReport analyze(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        int threads = 12,
        bool exclude_signature_123 = false,
        bool exclude_triple_family = false);
};

struct LocalSldProjectiveCoherenceCliOptions {
    std::string state_path;
    std::string certificate_path;
    int threads = 12;
    bool exclude_signature_123 = false;
    bool exclude_triple_family = false;
};

class LocalSldProjectiveCoherenceCli {
public:
    [[nodiscard]] static LocalSldProjectiveCoherenceCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const LocalSldProjectiveCoherenceCliOptions& options,
        std::ostream& out);
};

}  // namespace lemma
