// Consistency check against the existing C++ owners, not a Lean proof premise.
#include "local_sld_projective_normalization_objective.hpp"
#include "proof_scaling.hpp"
#include "local_sld_projective_normalization_cauchy_objective.hpp"
#include "local_sld_projective_normalization_schur_objective.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace lemma;

namespace {
void check(const char* name, SpectralReal actual, SpectralReal expected) {
    const SpectralReal tolerance = 1e-14L * std::max(std::abs(expected), 1e-20L);
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(std::string("C++/Lean witness mismatch: ") + name);
    }
}

SpectralState witness(int cutoff, SpectralReal amplitude) {
    SpectralState state;
    for (int x = -cutoff; x <= cutoff; ++x) {
        for (int y = -cutoff; y <= cutoff; ++y) {
            for (int z = -cutoff; z <= cutoff; ++z) {
                if (x == 0 && y == 0 && z == 0) continue;
                const WaveVector wave{x, y, z};
                state.index.emplace(wave, state.waves.size());
                state.waves.push_back(wave);
                state.velocity.push_back({});
            }
        }
    }
    const WaveVector waves[] = {{1, 0, 0}, {0, 1, 1}, {1, 1, 1},
                               {-1, 0, 0}, {0, -1, -1}, {-1, -1, -1}};
    const int imaginary[][3] = {{0, 1, 0}, {1, 0, 0}, {1, -1, 0},
                                {0, -1, 0}, {-1, 0, 0}, {-1, 1, 0}};
    for (std::size_t k = 0; k < 6; ++k) {
        auto& value = state.velocity.at(state.index.at(waves[k]));
        for (std::size_t j = 0; j < 3; ++j) {
            value[j] = SpectralComplex{0, amplitude * imaginary[k][j]};
        }
    }
    return state;
}
}  // namespace

int main() {
    try {
        SpectralGalerkin configuration;
        configuration.configure("direct", 1);
        const SpectralDynamics dynamics(configuration);
        const LocalSldProjectiveNormalizationObjective objective(
            dynamics, TriadPartition::local, 1);
        const LocalSldProjectiveNormalizationCauchyObjective cauchy(
            dynamics, TriadPartition::local, 1, 1);
        const LocalSldProjectiveNormalizationSchurObjective schur(
            dynamics, TriadPartition::local, 1, 2, 1);
        for (const int cutoff : {1, 2}) {
            for (const SpectralReal a : {1.0L, 0.5L, 0.25L}) {
                const SpectralState state = witness(cutoff, a);
                const auto result = objective.evaluate(state);
                const SpectralReal a2 = a * a;
                const SpectralReal a3 = a2 * a;
                check("E", SpectralStateOps::energy(state), 8 * a2);
                check("Z", result.enstrophy, 18 * a2);
                check("P", result.palinstrophy, 46 * a2);
                check("S_full", result.full_stretching, 2 * a3);
                check("S_selected", result.selected_stretching, 2 * a3);
                check("T_selected", result.selected_palinstrophy_cross, 6 * a3);
                check("S_core", result.fixed_core_stretching, 0);
                check("T_core", result.fixed_core_palinstrophy_cross, 0);
                check("PNT-2", result.palinstrophy_normalization_power_one,
                      1 / (876024 * a));
                if (cutoff == 1) {
                    const auto c = cauchy.evaluate(state);
                    const auto r = schur.evaluate(state);
                    const SpectralReal a4 = a2 * a2;
                    const SpectralReal majorant = 275.0L / (33366002112.0L * a2);
                    check("B1", c.selected_aggregate_h1_norm2, 22 * a4);
                    check("T2", c.tail_aggregate_h2_norm2, 50 * a4);
                    check("PNT-7", c.squared_cauchy_bound_power_one, majorant);
                    check("Schur B1", r.selected_aggregate_h1_norm2, 22 * a4);
                    check("D_H", r.diagonal_tail_h2_norm2, 50 * a4);
                    check("Gram row", r.normalized_absolute_gram_row_sum, 1);
                    check("PNT-12", r.schur_squared_majorant, majorant);
                    check("PNT-14", r.height_half_compensated_schur_squared_majorant, majorant);
                    std::cout << "C++ majorants PASS: B1=22a^4 T2=D=50a^4 R=1 PNT-7/8/12/14=275/(33366002112a^2)\n";
                }
                std::cout << "C++ witness PASS: cutoff=" << cutoff
                          << " amplitude=" << a << " PNT-2=1/(876024*amplitude)\n";
            }
        }
        const auto scaling = ScalingAnalyzer::analyze_palinstrophy_normalization();
        if (!(scaling.open_quotient_amplitude_degree == Rational(-1)) ||
            !(scaling.quotient_amplitude_degree == Rational(-2)) ||
            scaling.fixed_energy_search_covers_amplitude_ray) {
            throw std::runtime_error("C++/Lean amplitude-scope mismatch");
        }
        std::cout << "C++ amplitude/scope PASS: open=-1, squared=-2, fixed-energy excludes ray\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
