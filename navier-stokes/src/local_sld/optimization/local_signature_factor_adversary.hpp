#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSignatureFactorAdversaryOptions {
    int minimum_cutoff = 2;
    int maximum_cutoff = 6;
    int samples = 16;
    int workers = 12;
    SpectralReal viscosity = 0.1L;
    SpectralReal time_step = 0.0001L;
    std::uint64_t seed = 20260801;
    std::string certificate_path;
    std::string state_output_path;
};

struct LocalSignatureFactorAdversaryRow {
    int cutoff = 0;
    std::string profile;
    int valid_samples = 0;
    int opposite_direction_samples = 0;
    int simultaneous_growth_samples = 0;
    SpectralReal mean_log_derivative_product = 0.0L;
    SpectralReal maximum_log_derivative_product = 0.0L;
    SpectralReal maximum_simultaneous_growth_rate = 0.0L;
    SpectralReal maximum_critical_growth_rate = 0.0L;
    SpectralReal maximum_critical_density_derivative = 0.0L;
    SpectralReal maximum_simultaneous_initial_critical_density = 0.0L;
    SpectralReal maximum_simultaneous_critical_density_derivative = 0.0L;
    std::uint64_t product_seed = 0;
    std::uint64_t simultaneous_growth_seed = 0;
    std::uint64_t critical_growth_seed = 0;
    std::uint64_t critical_density_derivative_seed = 0;
    std::uint64_t simultaneous_density_seed = 0;
    std::uint64_t simultaneous_density_derivative_seed = 0;
};

struct LocalSignatureFactorAdversaryReport {
    int workers = 0;
    SpectralReal viscosity = 0.0L;
    SpectralReal time_step = 0.0L;
    bool universal_opposite_direction_rejected = false;
    bool universal_no_simultaneous_growth_rejected = false;
    std::vector<LocalSignatureFactorAdversaryRow> rows;
};

class LocalSignatureFactorAdversary {
public:
    [[nodiscard]] static LocalSignatureFactorAdversaryReport run(
        const LocalSignatureFactorAdversaryOptions& options);
};

class LocalSignatureFactorAdversaryCli {
public:
    [[nodiscard]] static LocalSignatureFactorAdversaryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const LocalSignatureFactorAdversaryOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
