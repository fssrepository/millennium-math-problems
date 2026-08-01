#include "state_analysis.hpp"

#include "spectral_galerkin.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace lemma {

SpectralState SpectralStateReader::read_tsv(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open spectral state: " + path);
    }
    std::string line;
    if (!std::getline(input, line) || !std::getline(input, line)) {
        throw std::runtime_error("spectral state header is incomplete: " + path);
    }
    SpectralState state;
    WaveVector wave;
    SpectralReal ux_re = 0.0L;
    SpectralReal ux_im = 0.0L;
    SpectralReal uy_re = 0.0L;
    SpectralReal uy_im = 0.0L;
    SpectralReal uz_re = 0.0L;
    SpectralReal uz_im = 0.0L;
    while (input >> wave.x >> wave.y >> wave.z >> ux_re >> ux_im >> uy_re >>
           uy_im >> uz_re >> uz_im) {
        if (norm_squared(wave) == 0) {
            throw std::runtime_error("spectral state contains the zero mode");
        }
        if (state.index.contains(wave)) {
            throw std::runtime_error("spectral state contains a duplicate wave");
        }
        state.index.emplace(wave, state.waves.size());
        state.waves.push_back(wave);
        state.velocity.push_back(
            {SpectralComplex{ux_re, ux_im}, SpectralComplex{uy_re, uy_im},
             SpectralComplex{uz_re, uz_im}});
    }
    if (state.waves.empty() || !input.eof()) {
        throw std::runtime_error("spectral state data is invalid: " + path);
    }
    for (const WaveVector mode : state.waves) {
        if (!state.index.contains(-mode)) {
            throw std::runtime_error(
                "spectral state is missing a conjugate wave");
        }
    }
    return state;
}

StateAnalysisReport SpectralStateAnalyzer::analyze(
    const SpectralState& state, const SpectralObjective& objective,
    const StateAnalysisOptions& options) {
    if (options.top_modes < 0 ||
        !(options.active_relative_tolerance >= 0.0L) ||
        !std::isfinite(options.active_relative_tolerance)) {
        throw std::invalid_argument("invalid state-analysis options");
    }
    StateAnalysisReport report;
    report.state_path = options.state_path;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.modes = static_cast<int>(state.waves.size());
    report.objective = objective.evaluate(state);
    report.active_threshold =
        options.active_relative_tolerance * report.objective.energy;
    report.shells.resize(static_cast<std::size_t>(report.cutoff + 1));
    for (int shell = 0; shell <= report.cutoff; ++shell) {
        report.shells[static_cast<std::size_t>(shell)].shell = shell;
    }
    const SpectralIncrement gradient =
        objective.energy_level_gradient(state);
    SpectralIncrement projected_gradient = gradient;
    SpectralReal radial_pairing = 0.0L;
    SpectralReal gradient_norm2 = 0.0L;
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        radial_pairing += std::real(dot_hermitian(
            projected_gradient[index], state.velocity[index]));
    }
    const SpectralReal radial_coefficient =
        radial_pairing / report.objective.energy;
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        for (std::size_t component = 0; component < 3; ++component) {
            projected_gradient[index][component] -=
                radial_coefficient * state.velocity[index][component];
        }
        gradient_norm2 += std::real(dot_hermitian(
            projected_gradient[index], projected_gradient[index]));
    }
    report.projected_q_gradient_norm =
        std::sqrt(std::max(0.0L, gradient_norm2));
    if (report.projected_q_gradient_norm > 0.0L) {
        for (ComplexVector& value : projected_gradient) {
            for (SpectralComplex& component : value) {
                component /= report.projected_q_gradient_norm;
            }
        }
        auto retracted_q = [&](SpectralReal step) {
            SpectralState candidate = state;
            for (std::size_t index = 0;
                 index < candidate.velocity.size(); ++index) {
                for (std::size_t component = 0; component < 3; ++component) {
                    candidate.velocity[index][component] +=
                        step * projected_gradient[index][component];
                }
            }
            SpectralStateOps::normalize_energy(
                candidate, report.objective.energy);
            return objective.evaluate(candidate).energy_level_quantity;
        };
        constexpr SpectralReal derivative_step = 1e-6L;
        report.retraction_directional_derivative =
            (retracted_q(derivative_step) -
             retracted_q(-derivative_step)) /
            (2.0L * derivative_step);
        report.retraction_gradient_relative_error =
            std::abs(report.retraction_directional_derivative -
                     report.projected_q_gradient_norm) /
            std::max(1e-30L, report.projected_q_gradient_norm);
        SpectralReal line_step = 0.2L;
        for (int sample = 0; sample < 16; ++sample) {
            report.q_line_profile.push_back(
                {line_step, retracted_q(line_step)});
            line_step *= 0.5L;
        }
    }
    std::vector<StateModeAnalysis> modes;
    modes.reserve(state.waves.size());
    for (std::size_t index = 0; index < state.waves.size(); ++index) {
        const WaveVector wave = state.waves[index];
        const int shell =
            std::max({std::abs(wave.x), std::abs(wave.y), std::abs(wave.z)});
        StateShellAnalysis& shell_report =
            report.shells[static_cast<std::size_t>(shell)];
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(wave));
        const SpectralReal mode_energy = std::real(
            dot_hermitian(state.velocity[index], state.velocity[index]));
        const SpectralReal gradient_norm2 =
            std::real(dot_hermitian(gradient[index], gradient[index]));
        ++shell_report.modes;
        shell_report.energy += mode_energy;
        shell_report.enstrophy += wave2 * mode_energy;
        shell_report.palinstrophy += wave2 * wave2 * mode_energy;
        shell_report.q_gradient_norm2 += gradient_norm2;
        if (mode_energy > report.active_threshold) {
            ++shell_report.active_modes;
            ++report.active_modes;
            report.highest_active_shell =
                std::max(report.highest_active_shell, shell);
        }
        const SpectralComplex divergence =
            wave_dot(wave, state.velocity[index]);
        const SpectralReal divergence_scale =
            std::sqrt(std::max(1e-30L, wave2 * mode_energy));
        report.divergence_residual = std::max(
            report.divergence_residual,
            std::abs(divergence) / divergence_scale);
        const std::size_t negative = state.index.at(-wave);
        SpectralReal reality_error2 = 0.0L;
        for (std::size_t component = 0; component < 3; ++component) {
            reality_error2 += std::norm(
                state.velocity[negative][component] -
                std::conj(state.velocity[index][component]));
        }
        report.reality_residual = std::max(
            report.reality_residual,
            std::sqrt(reality_error2 / std::max(1e-30L, mode_energy)));
        modes.push_back({wave, mode_energy, gradient_norm2});
    }
    std::sort(modes.begin(), modes.end(),
              [](const StateModeAnalysis& left,
                 const StateModeAnalysis& right) {
                  return left.energy > right.energy;
              });
    if (modes.size() > static_cast<std::size_t>(options.top_modes)) {
        modes.resize(static_cast<std::size_t>(options.top_modes));
    }
    report.top_modes = std::move(modes);
    return report;
}

void StateAnalysisReporter::write_console(const StateAnalysisReport& report,
                                          std::ostream& out) {
    out << std::setprecision(12)
        << "state=" << report.state_path << '\n'
        << "cutoff=" << report.cutoff << " modes=" << report.modes
        << " active_modes=" << report.active_modes
        << " highest_active_shell=" << report.highest_active_shell << '\n'
        << "E=" << static_cast<double>(report.objective.energy)
        << " Z=" << static_cast<double>(report.objective.enstrophy)
        << " P=" << static_cast<double>(report.objective.palinstrophy)
        << " Q="
        << static_cast<double>(report.objective.energy_level_quantity)
        << " divergence_residual="
        << static_cast<double>(report.divergence_residual)
        << " reality_residual=" << static_cast<double>(report.reality_residual)
        << " projected_q_gradient_norm="
        << static_cast<double>(report.projected_q_gradient_norm)
        << " retraction_gradient_error="
        << static_cast<double>(report.retraction_gradient_relative_error)
        << "\n\nshell,modes,active,energy,energy_fraction,enstrophy,palinstrophy,"
           "q_gradient_norm\n";
    for (const StateShellAnalysis& shell : report.shells) {
        out << shell.shell << ',' << shell.modes << ',' << shell.active_modes
            << ',' << static_cast<double>(shell.energy) << ','
            << static_cast<double>(
                   shell.energy / std::max(1e-30L, report.objective.energy))
            << ',' << static_cast<double>(shell.enstrophy) << ','
            << static_cast<double>(shell.palinstrophy) << ','
            << static_cast<double>(std::sqrt(shell.q_gradient_norm2)) << '\n';
    }
    out << "\ntop_mode,kx,ky,kz,energy,q_gradient_norm\n";
    for (std::size_t rank = 0; rank < report.top_modes.size(); ++rank) {
        const StateModeAnalysis& mode = report.top_modes[rank];
        out << rank + 1 << ',' << mode.wave.x << ',' << mode.wave.y << ','
            << mode.wave.z << ',' << static_cast<double>(mode.energy) << ','
            << static_cast<double>(std::sqrt(mode.q_gradient_norm2)) << '\n';
    }
    out << "\nline_step,Q,improvement\n";
    for (const StateLineSample& sample : report.q_line_profile) {
        out << static_cast<double>(sample.step) << ','
            << static_cast<double>(sample.q) << ','
            << static_cast<double>(sample.q -
                                   report.objective.energy_level_quantity)
            << '\n';
    }
}

void StateAnalysisReporter::write_json(const StateAnalysisReport& report,
                                       std::ostream& out) {
    out << std::setprecision(18)
        << "{\n  \"schema\": \"navier-stokes-state-analysis-v1\",\n"
        << "  \"state\": \"" << report.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff
        << ", \"modes\": " << report.modes
        << ", \"active_modes\": " << report.active_modes
        << ", \"highest_active_shell\": " << report.highest_active_shell
        << ",\n  \"E\": " << static_cast<double>(report.objective.energy)
        << ", \"Z\": " << static_cast<double>(report.objective.enstrophy)
        << ", \"P\": " << static_cast<double>(report.objective.palinstrophy)
        << ", \"Q\": "
        << static_cast<double>(report.objective.energy_level_quantity)
        << ",\n  \"divergence_residual\": "
        << static_cast<double>(report.divergence_residual)
        << ", \"reality_residual\": "
        << static_cast<double>(report.reality_residual)
        << ", \"projected_q_gradient_norm\": "
        << static_cast<double>(report.projected_q_gradient_norm)
        << ", \"retraction_directional_derivative\": "
        << static_cast<double>(report.retraction_directional_derivative)
        << ", \"retraction_gradient_relative_error\": "
        << static_cast<double>(report.retraction_gradient_relative_error)
        << ",\n  \"shells\": [\n";
    for (std::size_t index = 0; index < report.shells.size(); ++index) {
        const StateShellAnalysis& shell = report.shells[index];
        out << "    {\"shell\": " << shell.shell
            << ", \"modes\": " << shell.modes
            << ", \"active_modes\": " << shell.active_modes
            << ", \"energy\": " << static_cast<double>(shell.energy)
            << ", \"enstrophy\": " << static_cast<double>(shell.enstrophy)
            << ", \"palinstrophy\": "
            << static_cast<double>(shell.palinstrophy)
            << ", \"q_gradient_norm\": "
            << static_cast<double>(std::sqrt(shell.q_gradient_norm2)) << '}'
            << (index + 1 == report.shells.size() ? "\n" : ",\n");
    }
    out << "  ]\n}\n";
}

int run_state_analysis(const StateAnalysisOptions& options, std::ostream& out) {
    if (options.state_path.empty()) {
        throw std::invalid_argument("--state is required");
    }
    const SpectralState state = SpectralStateReader::read_tsv(options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("auto", 12);
    const SpectralDynamics dynamics(configuration);
    const SpectralObjective objective(dynamics);
    const StateAnalysisReport report =
        SpectralStateAnalyzer::analyze(state, objective, options);
    StateAnalysisReporter::write_console(report, out);
    if (!options.certificate_path.empty()) {
        std::ofstream certificate(options.certificate_path);
        if (!certificate) {
            throw std::runtime_error(
                "cannot open state-analysis certificate: " +
                options.certificate_path);
        }
        StateAnalysisReporter::write_json(report, certificate);
        out << "Certificate written to " << options.certificate_path << '\n';
    }
    return report.divergence_residual < 1e-12L &&
                   report.reality_residual < 1e-12L
               ? 0
               : 2;
}

StateAnalysisOptions StateAnalysisCli::parse(int argc, char** argv, int first) {
    StateAnalysisOptions options;
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
        } else if (name == "--certificate") {
            options.certificate_path = next_value(index, name);
        } else if (name == "--top") {
            options.top_modes = std::stoi(next_value(index, name));
        } else if (name == "--active-relative") {
            options.active_relative_tolerance =
                std::stold(next_value(index, name));
        } else {
            throw std::invalid_argument(
                "unknown state-analysis option: " + name);
        }
    }
    return options;
}

void StateAnalysisCli::print_help(std::ostream& out) {
    out << "State analysis options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write shell analysis JSON\n"
        << "  --top N               report N strongest modes (default 12)\n"
        << "  --active-relative X   relative energy activity threshold\n";
}

}  // namespace lemma
