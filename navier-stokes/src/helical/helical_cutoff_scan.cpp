#include "helical_cutoff_scan.hpp"

#include "helical_sector_adjoint.hpp"
#include "lemma_adversary.hpp"
#include "spectral_dynamics.hpp"
#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

HelicalSectorSelection selection_from_name(const std::string& name) {
    if (name == "heterochiral") {
        return HelicalSectorSelection::heterochiral();
    }
    if (name == "homochiral") {
        return HelicalSectorSelection::homochiral();
    }
    throw std::invalid_argument(
        "--selection must be heterochiral or homochiral");
}

HelicalLocalSpread spread_from_name(const std::string& name) {
    if (name == "all") {
        return HelicalLocalSpread::all;
    }
    if (name == "equal") {
        return HelicalLocalSpread::equal;
    }
    if (name == "narrow") {
        return HelicalLocalSpread::narrow;
    }
    if (name == "broad") {
        return HelicalLocalSpread::broad;
    }
    throw std::invalid_argument(
        "--spread must be all, equal, narrow, or broad");
}

void create_parent(const std::string& path) {
    const std::filesystem::path parent =
        std::filesystem::path(path).parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

SpectralReal relative_difference(SpectralReal left, SpectralReal right) {
    const SpectralReal scale = std::max(
        {std::abs(left), std::abs(right),
         std::numeric_limits<SpectralReal>::min()});
    return std::abs(left - right) / scale;
}

SpectralReal cancellation_ratio(
    SpectralReal signed_value, SpectralReal absolute_sum) {
    return absolute_sum > 0.0L
        ? std::abs(signed_value) / absolute_sum
        : 0.0L;
}

}  // namespace

HelicalCutoffScanOptions HelicalCutoffScan::parse(
    int argc, char** argv, int first) {
    HelicalCutoffScanOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_path = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--spread") {
            options.spread = next(index, name);
        } else if (name == "--min-cutoff") {
            options.minimum_cutoff = std::stoi(next(index, name));
        } else if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next(index, name));
        } else if (name == "--workers") {
            options.workers = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else if (name == "--tolerance") {
            options.convergence_tolerance = std::stold(next(index, name));
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown helical-cutoff-scan option: " + name);
        }
    }
    return options;
}

void HelicalCutoffScan::print_help(std::ostream& out) {
    out << "Helical cutoff-scan options:\n"
        << "  --state PATH          input replayable Fourier-state TSV\n"
        << "  --certificate PATH    write aggregate cutoff JSON\n"
        << "  --selection NAME      heterochiral or homochiral\n"
        << "  --spread NAME         all, equal, narrow, or broad local triads\n"
        << "  --min-cutoff K        first Galerkin cutoff\n"
        << "  --max-cutoff K        last Galerkin cutoff\n"
        << "  --trajectory-steps N  RK4 steps at every cutoff\n"
        << "  --workers N           parallel cutoff workers\n"
        << "  --nu X                viscosity\n"
        << "  --dt X                RK4 time step\n"
        << "  --tolerance X         last-pair relative tolerance\n"
        << "  --seed N              deterministic lift seed\n";
}

int HelicalCutoffScan::run(
    const HelicalCutoffScanOptions& options, std::ostream& out) {
    if (options.state_path.empty() || options.certificate_path.empty()) {
        throw std::invalid_argument(
            "helical-cutoff-scan requires --state and --certificate");
    }
    if (options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.trajectory_steps < 1 || options.workers < 1 ||
        !(options.viscosity > 0.0L) || !(options.time_step > 0.0L) ||
        !(options.convergence_tolerance > 0.0L)) {
        throw std::invalid_argument("invalid helical cutoff-scan options");
    }
    const HelicalSectorSelection selection =
        selection_from_name(options.selection).with_spread(
            spread_from_name(options.spread));
    const SpectralState base = SpectralStateReader::read_tsv(
        options.state_path);
    const int base_cutoff = SpectralStateOps::cutoff(base);
    const SpectralReal base_energy = SpectralStateOps::energy(base);
    const std::size_t count = static_cast<std::size_t>(
        options.maximum_cutoff - options.minimum_cutoff + 1);
    std::vector<SpectralState> states;
    states.reserve(count);
    for (int cutoff = options.minimum_cutoff;
         cutoff <= options.maximum_cutoff; ++cutoff) {
        if (cutoff < base_cutoff) {
            states.push_back(SpectralStateFactory::project(base, cutoff));
        } else if (cutoff == base_cutoff) {
            states.push_back(base);
        } else {
            std::mt19937_64 generator(
                options.seed + static_cast<std::uint64_t>(cutoff) *
                    UINT64_C(0x9e3779b97f4a7c15));
            states.push_back(SpectralStateFactory::lift(
                base, cutoff, generator));
        }
        SpectralStateOps::normalize_energy(states.back(), base_energy);
        static_cast<void>(SpectralStateOps::interactions(states.back()));
    }

    std::vector<SpectralReal> objectives(count, 0.0L);
    std::vector<HelicalTriadReport> initial_helicity(count);
    std::vector<HelicalTriadReport> final_helicity(count);
    const LemmaAdversary executor(options.workers);
    executor.run_restarts(count, [&](std::size_t index) {
        SpectralGalerkin galerkin;
        galerkin.configure("direct", 1);
        const SpectralDynamics dynamics(galerkin);
        const HelicalSectorAdjoint adjoint(dynamics);
        initial_helicity[index] = HelicalTriadLedger::analyze(states[index]);
        const HelicalSectorTrajectoryValue trajectory =
            adjoint.critical_trajectory(
                states[index], options.viscosity, options.time_step,
                options.trajectory_steps, selection);
        objectives[index] = trajectory.objective_value;
        final_helicity[index] = HelicalTriadLedger::analyze(
            trajectory.final_state);
    });

    std::vector<SpectralReal> adjacent_relative(count, 0.0L);
    for (std::size_t index = 1; index < count; ++index) {
        adjacent_relative[index] = relative_difference(
            objectives[index - 1], objectives[index]);
    }
    const SpectralReal last_relative = count > 1
        ? adjacent_relative.back()
        : 0.0L;
    const bool converged = count > 1 &&
        last_relative <= options.convergence_tolerance;
    const SpectralReal horizon =
        options.time_step * static_cast<SpectralReal>(options.trajectory_steps);

    create_parent(options.certificate_path);
    std::ofstream certificate(options.certificate_path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write cutoff certificate: " +
            options.certificate_path);
    }
    certificate << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-helical-cutoff-scan-v1\",\n"
        << "  \"input_state\": \"" << options.state_path << "\",\n"
        << "  \"selection\": \"" << options.selection << "\",\n"
        << "  \"spread\": \"" << options.spread << "\",\n"
        << "  \"minimum_cutoff\": " << options.minimum_cutoff << ",\n"
        << "  \"maximum_cutoff\": " << options.maximum_cutoff << ",\n"
        << "  \"workers\": " << executor.threads() << ",\n"
        << "  \"trajectory_steps\": " << options.trajectory_steps << ",\n"
        << "  \"viscosity\": " << static_cast<double>(options.viscosity)
        << ",\n  \"time_step\": " << static_cast<double>(options.time_step)
        << ",\n  \"time_horizon\": " << static_cast<double>(horizon)
        << ",\n  \"convergence_tolerance\": "
        << static_cast<double>(options.convergence_tolerance)
        << ",\n  \"last_adjacent_relative_difference\": "
        << static_cast<double>(last_relative)
        << ",\n  \"numerically_converged\": "
        << (converged ? "true" : "false") << ",\n"
        << "  \"samples\": [\n";
    for (std::size_t index = 0; index < count; ++index) {
        const HelicalTriadReport& initial = initial_helicity[index];
        const HelicalTriadReport& final = final_helicity[index];
        certificate << "    {\"cutoff\": "
            << options.minimum_cutoff + static_cast<int>(index)
            << ", \"objective\": " << static_cast<double>(objectives[index])
            << ", \"mean_density\": "
            << static_cast<double>(objectives[index] / horizon)
            << ", \"previous_relative_difference\": "
            << static_cast<double>(adjacent_relative[index])
            << ", \"initial_positive_helical_energy\": "
            << static_cast<double>(initial.positive_helical_energy)
            << ", \"initial_negative_helical_energy\": "
            << static_cast<double>(initial.negative_helical_energy)
            << ", \"initial_heterochiral_signed_local\": "
            << static_cast<double>(initial.heterochiral_local_stretching)
            << ", \"initial_heterochiral_absolute_local\": "
            << static_cast<double>(
                   initial.heterochiral_absolute_local_stretching)
            << ", \"initial_heterochiral_cancellation_ratio\": "
            << static_cast<double>(cancellation_ratio(
                   initial.heterochiral_local_stretching,
                   initial.heterochiral_absolute_local_stretching))
            << ", \"final_positive_helical_energy\": "
            << static_cast<double>(final.positive_helical_energy)
            << ", \"final_negative_helical_energy\": "
            << static_cast<double>(final.negative_helical_energy)
            << ", \"final_heterochiral_signed_local\": "
            << static_cast<double>(final.heterochiral_local_stretching)
            << ", \"final_heterochiral_absolute_local\": "
            << static_cast<double>(
                   final.heterochiral_absolute_local_stretching)
            << ", \"final_heterochiral_cancellation_ratio\": "
            << static_cast<double>(cancellation_ratio(
                   final.heterochiral_local_stretching,
                   final.heterochiral_absolute_local_stretching)) << "}"
            << (index + 1 == count ? "\n" : ",\n");
    }
    certificate << "  ]\n}\n";

    out << std::setprecision(12)
        << "helical cutoff scan " << options.selection
        << " spread=" << options.spread
        << " K=" << options.minimum_cutoff << ".."
        << options.maximum_cutoff << " horizon="
        << static_cast<double>(horizon) << '\n';
    for (std::size_t index = 0; index < count; ++index) {
        out << "  K=" << options.minimum_cutoff + static_cast<int>(index)
            << " objective=" << static_cast<double>(objectives[index])
            << " relative="
            << static_cast<double>(adjacent_relative[index]) << '\n';
    }
    out << "last-pair convergence=" << (converged ? "PASS" : "FAIL")
        << " tolerance="
        << static_cast<double>(options.convergence_tolerance) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return converged ? 0 : 2;
}

}  // namespace lemma
