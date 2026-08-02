#pragma once

#include "spectral_state.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct LocalSignatureAdversaryOptions {
    int minimum_cutoff = 1;
    int maximum_cutoff = 5;
    int samples = 24;
    int workers = 12;
    std::uint64_t seed = 20260801;
    std::string certificate_path;
};

struct LocalSignatureAdversaryRow {
    int cutoff = 0;
    std::string profile;
    int samples = 0;
    SpectralReal mean_effective_count = 0.0L;
    SpectralReal maximum_effective_count = 0.0L;
    SpectralReal maximizer_dominant_fraction = 0.0L;
    std::size_t maximum_coherent_signatures = 0;
    std::uint64_t maximizing_seed = 0;
    SpectralReal mean_signed_amplification = 0.0L;
    SpectralReal maximum_signed_amplification = 0.0L;
    std::uint64_t amplifying_seed = 0;
};

struct LocalSignatureAdversaryReport {
    int workers = 0;
    SpectralReal maximum_observed_effective_count = 0.0L;
    SpectralReal flat_maximum_log_slope = 0.0L;
    SpectralReal outer_maximum_log_slope = 0.0L;
    SpectralReal flat_amplification_log_slope = 0.0L;
    SpectralReal outer_amplification_log_slope = 0.0L;
    bool finite_state_search_is_not_a_proof = true;
    std::vector<LocalSignatureAdversaryRow> rows;
};

class LocalSignatureAdversary {
public:
    [[nodiscard]] static LocalSignatureAdversaryReport run(
        const LocalSignatureAdversaryOptions& options);
};

class LocalSignatureAdversaryCli {
public:
    [[nodiscard]] static LocalSignatureAdversaryOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const LocalSignatureAdversaryOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
