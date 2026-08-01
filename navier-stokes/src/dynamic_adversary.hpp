#pragma once

#include "gradient_adversary.hpp"
#include "trajectory_analyzer.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace lemma {

struct DynamicAdversaryOptions {
    int generations = 0;
    SpectralReal mutation = 0.2L;
    SpectralReal viscosity = 0.1L;
    SpectralReal final_time = 0.1L;
    SpectralReal time_step = 0.002L;
    std::uint64_t seed = 20260801;
    std::string objective = "critical-integral";
    std::string optimizer = "gradient";
    std::string gradient_method = "steepest";
    int sobolev_order = 0;
    SpectralReal sobolev_cap = 0.0L;
    int minimum_dyadic_gap = 2;
};

struct DynamicAdversaryResult {
    SpectralState state;
    StaticObjective initial_objective;
    EvolutionResult evolution;
    EvolutionResult refined_evolution;
    SpectralReal time_step_relative_error = 0.0L;
    SpectralReal search_initial_objective = 0.0L;
    SpectralReal search_final_objective = 0.0L;
    std::vector<GradientIterationRecord> gradient_trace;
    std::vector<SpectralReal> restart_objectives;
    int winning_restart = 0;
    int accepted_mutations = 0;
    int accepted_gradient_steps = 0;
    int evaluations = 0;
};

class DynamicAdversary {
public:
    DynamicAdversary(std::string backend, int compute_threads);

    [[nodiscard]] DynamicAdversaryResult optimize(
        const SpectralState& primary_start,
        const SpectralState* secondary_start,
        const DynamicAdversaryOptions& options,
        bool refine = true) const;

    void refine(DynamicAdversaryResult& result,
                const DynamicAdversaryOptions& options) const;

    [[nodiscard]] static SpectralReal objective_value(
        const EvolutionResult& evolution, const std::string& objective);

private:
    SpectralGalerkin galerkin_;
    SpectralDynamics dynamics_;
    SpectralObjective objective_;
    TrajectoryAnalyzer trajectory_;
    SpectralAdjoint adjoint_;
    GradientAdversary gradient_;
};

class DynamicAdversaryEnsemble {
public:
    DynamicAdversaryEnsemble(std::string backend, int workers);

    [[nodiscard]] DynamicAdversaryResult optimize(
        const SpectralState& primary_start,
        const SpectralState* secondary_start,
        const DynamicAdversaryOptions& options,
        int restarts) const;

private:
    std::string backend_;
    int workers_ = 1;
};

}  // namespace lemma
