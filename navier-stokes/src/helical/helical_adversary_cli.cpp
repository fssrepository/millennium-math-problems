#include "helical_adversary_cli.hpp"

#include "helical_sector_adjoint.hpp"
#include "helical_sector_adversary.hpp"
#include "helical_trajectory_adversary.hpp"
#include "lemma_adversary.hpp"
#include "spectral_dynamics.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <filesystem>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

HelicalSectorSelection parse_selection(const std::string& selection) {
    if (selection == "heterochiral") {
        return HelicalSectorSelection::heterochiral();
    }
    if (selection == "homochiral") {
        return HelicalSectorSelection::homochiral();
    }
    throw std::invalid_argument(
        "--selection must be homochiral or heterochiral");
}

HelicalLocalSpread parse_spread(const std::string& spread) {
    if (spread == "all") {
        return HelicalLocalSpread::all;
    }
    if (spread == "equal") {
        return HelicalLocalSpread::equal;
    }
    if (spread == "narrow") {
        return HelicalLocalSpread::narrow;
    }
    if (spread == "broad") {
        return HelicalLocalSpread::broad;
    }
    throw std::invalid_argument(
        "--spread must be all, equal, narrow, or broad");
}

void ensure_parent_directory(const std::string& path) {
    const std::filesystem::path parent =
        std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

void write_state(const SpectralState& state, const std::string& path,
                 const std::string& mode, const std::string& selection,
                 const std::string& spread, SpectralReal objective) {
    ensure_parent_directory(path);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot write helical state: " + path);
    }
    output << std::setprecision(21)
           << "# cutoff=" << SpectralStateOps::cutoff(state)
           << " energy=" << static_cast<double>(
                  SpectralStateOps::energy(state))
           << " mode=" << mode << " selection=" << selection
           << " spread=" << spread
           << " objective=" << static_cast<double>(objective) << '\n'
           << "kx\tky\tkz\tux_re\tux_im\tuy_re\tuy_im\tuz_re\tuz_im\n";
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        output << wave.x << '\t' << wave.y << '\t' << wave.z;
        for (const SpectralComplex component : state.velocity[index]) {
            output << '\t' << static_cast<double>(component.real())
                   << '\t' << static_cast<double>(component.imag());
        }
        output << '\n';
    }
}

}  // namespace

HelicalAdversaryCliOptions HelicalAdversaryCli::parse(
    int argc, char** argv, int first) {
    HelicalAdversaryCliOptions options;
    auto next_value = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_path = next_value(index, name);
        } else if (name == "--state-output") {
            options.state_output_path = next_value(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next_value(index, name);
        } else if (name == "--selection") {
            options.selection = next_value(index, name);
        } else if (name == "--mode") {
            options.mode = next_value(index, name);
        } else if (name == "--spread") {
            options.spread = next_value(index, name);
        } else if (name == "--cutoff") {
            options.cutoff = std::stoi(next_value(index, name));
        } else if (name == "--iterations") {
            options.iterations = std::stoi(next_value(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next_value(index, name));
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next_value(index, name));
        } else if (name == "--restarts") {
            options.restarts = std::stoi(next_value(index, name));
        } else if (name == "--workers") {
            options.workers = std::stoi(next_value(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next_value(index, name));
        } else if (name == "--restart-mutation") {
            options.restart_mutation = std::stold(next_value(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next_value(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next_value(index, name));
        } else if (name == "--seed") {
            options.seed = std::stoull(next_value(index, name));
        } else {
            throw std::invalid_argument(
                "unknown helical-adversary option: " + name);
        }
    }
    return options;
}

void HelicalAdversaryCli::print_help(std::ostream& out) {
    out << "Helical adversary options:\n"
        << "  --state PATH          input replayable Fourier-state TSV\n"
        << "  --state-output PATH   write optimized Fourier-state TSV\n"
        << "  --certificate PATH    write optimization JSON\n"
        << "  --selection NAME      heterochiral or homochiral\n"
        << "  --mode NAME           trajectory or static\n"
        << "  --spread NAME         all, equal, narrow, or broad local triads\n"
        << "  --cutoff K            project/lift input to cutoff K\n"
        << "  --iterations N        exact-gradient iterations (0 evaluates only)\n"
        << "  --line-search N       line-search trials per iteration\n"
        << "  --trajectory-steps N  RK4 steps for trajectory mode\n"
        << "  --restarts N          independent starts (first is input)\n"
        << "  --workers N           parallel restart workers\n"
        << "  --step X              initial Riemannian step\n"
        << "  --restart-mutation X  perturbation size for extra starts\n"
        << "  --nu X                viscosity\n"
        << "  --dt X                RK4 time step\n"
        << "  --seed N              deterministic lift seed\n";
}

int run_helical_adversary(
    const HelicalAdversaryCliOptions& options, std::ostream& out) {
    if (options.state_path.empty() || options.state_output_path.empty() ||
        options.certificate_path.empty()) {
        throw std::invalid_argument(
            "helical-adversary requires --state, --state-output, and --certificate");
    }
    if (options.mode != "static" && options.mode != "trajectory") {
        throw std::invalid_argument("--mode must be static or trajectory");
    }
    if (options.iterations < 0 || options.restarts < 1 || options.workers < 1 ||
        options.restart_mutation < 0.0L) {
        throw std::invalid_argument(
            "iterations/mutation must be nonnegative; restarts/workers positive");
    }
    const HelicalSectorSelection selection =
        parse_selection(options.selection).with_spread(
            parse_spread(options.spread));
    SpectralState state = SpectralStateReader::read_tsv(options.state_path);
    const SpectralReal energy = SpectralStateOps::energy(state);
    const int input_cutoff = SpectralStateOps::cutoff(state);
    const int target_cutoff = options.cutoff == 0
        ? input_cutoff
        : options.cutoff;
    std::mt19937_64 generator(options.seed);
    if (target_cutoff < input_cutoff) {
        state = SpectralStateFactory::project(state, target_cutoff);
    } else if (target_cutoff > input_cutoff) {
        state = SpectralStateFactory::lift(state, target_cutoff, generator);
    }
    SpectralStateOps::normalize_energy(state, energy);

    std::vector<SpectralState> starts;
    starts.reserve(static_cast<std::size_t>(options.restarts));
    starts.push_back(state);
    for (int restart = 1; restart < options.restarts; ++restart) {
        std::mt19937_64 restart_generator(
            options.seed + static_cast<std::uint64_t>(restart) *
                UINT64_C(0x9e3779b97f4a7c15));
        starts.push_back(SpectralStateFactory::mutate(
            state, options.restart_mutation, restart_generator, false));
        SpectralStateOps::normalize_energy(starts.back(), energy);
    }
    static_cast<void>(SpectralStateOps::interactions(starts.front()));

    struct SearchResult {
        SpectralState state;
        SpectralReal initial_objective = 0.0L;
        SpectralReal objective = 0.0L;
        int accepted_steps = 0;
        int evaluations = 0;
    };
    std::vector<SearchResult> results(
        static_cast<std::size_t>(options.restarts));
    const LemmaAdversary executor(options.workers);
    SpectralState optimized;
    SpectralReal initial_objective = 0.0L;
    SpectralReal final_objective = 0.0L;
    int accepted_steps = 0;
    int evaluations = 0;
    if (options.iterations == 0) {
        executor.run_restarts(results.size(), [&](std::size_t restart) {
            SpectralReal value = 0.0L;
            if (options.mode == "static") {
                value = HelicalSectorObjective::evaluate(
                    starts[restart], selection).critical_integrand;
            } else {
                SpectralGalerkin local_galerkin;
                local_galerkin.configure("direct", 1);
                const SpectralDynamics local_dynamics(local_galerkin);
                const HelicalSectorAdjoint local_adjoint(local_dynamics);
                value = local_adjoint.critical_integral(
                    starts[restart], options.viscosity, options.time_step,
                    options.trajectory_steps, selection);
            }
            results[restart] = {starts[restart], value, value, 0, 1};
        });
    } else if (options.mode == "static") {
        HelicalSectorAdversaryOptions search;
        search.selection = selection;
        search.iterations = options.iterations;
        search.line_search_steps = options.line_search_steps;
        search.initial_step = options.initial_step;
        executor.run_restarts(results.size(), [&](std::size_t restart) {
            const HelicalSectorAdversaryResult result =
                HelicalSectorAdversary::maximize(starts[restart], search);
            results[restart] = {result.state, result.initial_objective,
                                result.objective, result.accepted_steps,
                                result.evaluations};
        });
    } else {
        HelicalTrajectoryAdversaryOptions search;
        search.selection = selection;
        search.iterations = options.iterations;
        search.line_search_steps = options.line_search_steps;
        search.trajectory_steps = options.trajectory_steps;
        search.initial_step = options.initial_step;
        search.viscosity = options.viscosity;
        search.time_step = options.time_step;
        executor.run_restarts(results.size(), [&](std::size_t restart) {
            SpectralGalerkin local_galerkin;
            local_galerkin.configure("direct", 1);
            const SpectralDynamics local_dynamics(local_galerkin);
            const HelicalSectorAdjoint local_adjoint(local_dynamics);
            const HelicalTrajectoryAdversaryResult result =
                HelicalTrajectoryAdversary::maximize(
                    starts[restart], search, local_adjoint);
            results[restart] = {result.state, result.initial_objective,
                                result.objective, result.accepted_steps,
                                result.evaluations};
        });
    }
    std::size_t winner = 0;
    int total_evaluations = 0;
    for (std::size_t restart = 0; restart < results.size(); ++restart) {
        total_evaluations += results[restart].evaluations;
        if (results[restart].objective > results[winner].objective) {
            winner = restart;
        }
    }
    optimized = std::move(results[winner].state);
    initial_objective = results[winner].initial_objective;
    final_objective = results[winner].objective;
    accepted_steps = results[winner].accepted_steps;
    evaluations = results[winner].evaluations;

    write_state(optimized, options.state_output_path, options.mode,
                options.selection, options.spread, final_objective);
    ensure_parent_directory(options.certificate_path);
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write helical certificate: " +
            options.certificate_path);
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-helical-adversary-v1\",\n"
        << "  \"input_state\": \"" << options.state_path << "\",\n"
        << "  \"output_state\": \"" << options.state_output_path << "\",\n"
        << "  \"mode\": \"" << options.mode << "\",\n"
        << "  \"selection\": \"" << options.selection << "\",\n"
        << "  \"spread\": \"" << options.spread << "\",\n"
        << "  \"cutoff\": " << target_cutoff << ",\n"
        << "  \"iterations\": " << options.iterations << ",\n"
        << "  \"restarts\": " << options.restarts << ",\n"
        << "  \"workers\": " << executor.threads() << ",\n"
        << "  \"winner_restart\": " << winner << ",\n"
        << "  \"restart_mutation\": "
        << static_cast<double>(options.restart_mutation) << ",\n"
        << "  \"accepted_steps\": " << accepted_steps << ",\n"
        << "  \"evaluations\": " << evaluations << ",\n"
        << "  \"total_evaluations\": " << total_evaluations << ",\n"
        << "  \"trajectory_steps\": "
        << (options.mode == "trajectory" ? options.trajectory_steps : 0)
        << ",\n  \"viscosity\": " << static_cast<double>(options.viscosity)
        << ",\n  \"time_step\": " << static_cast<double>(options.time_step)
        << ",\n  \"initial_objective\": "
        << static_cast<double>(initial_objective)
        << ",\n  \"final_objective\": "
        << static_cast<double>(final_objective)
        << ",\n  \"energy_error\": "
        << static_cast<double>(
               std::abs(SpectralStateOps::energy(optimized) - energy))
        << "\n}\n";
    out << std::setprecision(12)
        << "helical " << options.mode << " " << options.selection
        << " spread=" << options.spread
        << " cutoff=" << target_cutoff
        << " objective=" << static_cast<double>(initial_objective)
        << " -> " << static_cast<double>(final_objective)
        << " accepted=" << accepted_steps
        << " evaluations=" << evaluations
        << " restarts=" << options.restarts
        << " workers=" << executor.threads()
        << " winner=" << winner
        << " total_evaluations=" << total_evaluations << '\n'
        << "State written to " << options.state_output_path << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return options.iterations == 0 || final_objective > initial_objective
        ? 0
        : 2;
}

}  // namespace lemma
