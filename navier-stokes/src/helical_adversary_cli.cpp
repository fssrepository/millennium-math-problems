#include "helical_adversary_cli.hpp"

#include "helical_sector_adjoint.hpp"
#include "helical_sector_adversary.hpp"
#include "helical_trajectory_adversary.hpp"
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

void ensure_parent_directory(const std::string& path) {
    const std::filesystem::path parent =
        std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

void write_state(const SpectralState& state, const std::string& path,
                 const std::string& mode, const std::string& selection,
                 SpectralReal objective) {
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
        } else if (name == "--cutoff") {
            options.cutoff = std::stoi(next_value(index, name));
        } else if (name == "--iterations") {
            options.iterations = std::stoi(next_value(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next_value(index, name));
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next_value(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next_value(index, name));
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
        << "  --cutoff K            project/lift input to cutoff K\n"
        << "  --iterations N        exact-gradient iterations\n"
        << "  --line-search N       line-search trials per iteration\n"
        << "  --trajectory-steps N  RK4 steps for trajectory mode\n"
        << "  --step X              initial Riemannian step\n"
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
    const HelicalSectorSelection selection =
        parse_selection(options.selection);
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

    SpectralGalerkin galerkin;
    galerkin.configure("direct", 1);
    SpectralDynamics dynamics(galerkin);
    const HelicalSectorAdjoint adjoint(dynamics);
    SpectralState optimized;
    SpectralReal initial_objective = 0.0L;
    SpectralReal final_objective = 0.0L;
    int accepted_steps = 0;
    int evaluations = 0;
    if (options.mode == "static") {
        HelicalSectorAdversaryOptions search;
        search.selection = selection;
        search.iterations = options.iterations;
        search.line_search_steps = options.line_search_steps;
        search.initial_step = options.initial_step;
        const HelicalSectorAdversaryResult result =
            HelicalSectorAdversary::maximize(state, search);
        optimized = result.state;
        initial_objective = result.initial_objective;
        final_objective = result.objective;
        accepted_steps = result.accepted_steps;
        evaluations = result.evaluations;
    } else {
        HelicalTrajectoryAdversaryOptions search;
        search.selection = selection;
        search.iterations = options.iterations;
        search.line_search_steps = options.line_search_steps;
        search.trajectory_steps = options.trajectory_steps;
        search.initial_step = options.initial_step;
        search.viscosity = options.viscosity;
        search.time_step = options.time_step;
        const HelicalTrajectoryAdversaryResult result =
            HelicalTrajectoryAdversary::maximize(state, search, adjoint);
        optimized = result.state;
        initial_objective = result.initial_objective;
        final_objective = result.objective;
        accepted_steps = result.accepted_steps;
        evaluations = result.evaluations;
    }

    write_state(optimized, options.state_output_path, options.mode,
                options.selection, final_objective);
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
        << "  \"cutoff\": " << target_cutoff << ",\n"
        << "  \"iterations\": " << options.iterations << ",\n"
        << "  \"accepted_steps\": " << accepted_steps << ",\n"
        << "  \"evaluations\": " << evaluations << ",\n"
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
        << " cutoff=" << target_cutoff
        << " objective=" << static_cast<double>(initial_objective)
        << " -> " << static_cast<double>(final_objective)
        << " accepted=" << accepted_steps
        << " evaluations=" << evaluations << '\n'
        << "State written to " << options.state_output_path << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return final_objective > initial_objective ? 0 : 2;
}

}  // namespace lemma
