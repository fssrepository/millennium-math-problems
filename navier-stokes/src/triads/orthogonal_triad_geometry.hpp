#pragma once

#include "proof_scaling.hpp"
#include "spectral_state.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct OrthogonalTriadCountRow {
    int cutoff = 0;
    std::size_t ordered_pairs = 0;
    std::size_t maximum_input_degree = 0;
    std::size_t maximum_target_degree = 0;
    std::size_t elementary_degree_bound = 0;
};

struct OrthogonalTriadGeometryCertificate {
    int maximum_cutoff = 0;
    SpectralReal maximum_input_degree_ratio = 0.0L;
    SpectralReal maximum_target_degree_ratio = 0.0L;
    bool all_degree_bounds_hold = true;
    std::vector<OrthogonalTriadCountRow> rows;
};

struct OrthogonalTriadClosure {
    Rational derivative_weight_power{3};
    Rational generic_local_degree_power{3};
    Rational critical_degree_power{2};
    Rational interaction_degree_power{1};
    Rational cauchy_degree_power{1, 2};
    Rational transfer_frequency_power{7, 2};
    Rational generic_local_transfer_frequency_power{9, 2};
    Rational critical_transfer_frequency_power{4};
    Rational transfer_shell_energy_power{3, 2};
    Rational viscous_frequency_power{4};
    Rational viscous_shell_energy_power{1};
    Rational transfer_to_viscosity_frequency_power{-1, 2};
    Rational transfer_to_viscosity_energy_power{1, 2};
    bool high_frequency_absorbable_from_energy = true;
    bool orthogonal_degree_is_subcritical = true;
    bool generic_local_degree_is_supercritical = true;
};

class OrthogonalTriadGeometry {
public:
    [[nodiscard]] static OrthogonalTriadGeometryCertificate certify(
        int maximum_cutoff);
    [[nodiscard]] static OrthogonalTriadClosure analyze_closure();
};

struct OrthogonalTriadCliOptions {
    int maximum_cutoff = 8;
    std::string certificate_path;
};

class OrthogonalTriadCli {
public:
    [[nodiscard]] static OrthogonalTriadCliOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(const OrthogonalTriadCliOptions& options,
                   std::ostream& out);
};

}  // namespace lemma
