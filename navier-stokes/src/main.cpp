#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "lemma_engine.hpp"
#include "helical_adversary_cli.hpp"
#include "helical_cutoff_scan.hpp"
#include "state_analysis.hpp"
#include "orthogonal_triad_geometry.hpp"
#include "local_signature_geometry.hpp"
#include "local_signature_adversary.hpp"
#include "local_signature_factor_adversary.hpp"
#include "local_signature_gradient_adversary.hpp"
#include "local_signature_trajectory.hpp"
#include "local_quartic_closure_cli.hpp"
#include "local_sld_cyclic_ansatz.hpp"
#include "local_sld_cyclic_krylov_ansatz.hpp"
#include "local_sld_cyclic_trajectory_ansatz.hpp"
#include "local_sld_response_hierarchy.hpp"
#include "local_sld_signature_block.hpp"
#include "shifted_critical_density_cli.hpp"

namespace ns {

constexpr double pi = 3.141592653589793238462643383279502884;

struct Grid {
    int n;
    double length;
    double h;
    std::size_t cells;

    explicit Grid(int n_in, double length_in = 2.0 * pi)
        : n(n_in), length(length_in), h(length_in / static_cast<double>(n_in)),
          cells(static_cast<std::size_t>(n_in) * static_cast<std::size_t>(n_in) *
                static_cast<std::size_t>(n_in)) {
        if (n < 4) {
            throw std::invalid_argument("the grid size must be at least 4");
        }
    }

    [[nodiscard]] int wrap(int i) const {
        if (i < 0) {
            return i + n;
        }
        if (i >= n) {
            return i - n;
        }
        return i;
    }

    [[nodiscard]] std::size_t at(int i, int j, int k) const {
        const auto ii = static_cast<std::size_t>(wrap(i));
        const auto jj = static_cast<std::size_t>(wrap(j));
        const auto kk = static_cast<std::size_t>(wrap(k));
        const auto nn = static_cast<std::size_t>(n);
        return (kk * nn + jj) * nn + ii;
    }
};

struct Field {
    std::array<std::vector<double>, 3> c;

    explicit Field(std::size_t cells = 0) {
        for (auto& component : c) {
            component.assign(cells, 0.0);
        }
    }
};

struct ScalarField {
    std::vector<double> v;

    explicit ScalarField(std::size_t cells = 0) : v(cells, 0.0) {}
};

[[nodiscard]] double mean(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) /
           static_cast<double>(values.size());
}

void remove_mean(std::vector<double>& values) {
    const double m = mean(values);
    for (double& x : values) {
        x -= m;
    }
}

[[nodiscard]] double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t q = 0; q < a.size(); ++q) {
        result += a[q] * b[q];
    }
    return result;
}

void apply_negative_laplacian(const Grid& grid, const std::vector<double>& x,
                              std::vector<double>& result) {
    const double inv_h2 = 1.0 / (grid.h * grid.h);
    for (int k = 0; k < grid.n; ++k) {
        for (int j = 0; j < grid.n; ++j) {
            for (int i = 0; i < grid.n; ++i) {
                const auto q = grid.at(i, j, k);
                result[q] = (6.0 * x[q] - x[grid.at(i + 1, j, k)] -
                             x[grid.at(i - 1, j, k)] - x[grid.at(i, j + 1, k)] -
                             x[grid.at(i, j - 1, k)] - x[grid.at(i, j, k + 1)] -
                             x[grid.at(i, j, k - 1)]) *
                            inv_h2;
            }
        }
    }
}

struct CgResult {
    int iterations = 0;
    double relative_residual = 0.0;
    bool converged = true;
};

CgResult solve_periodic_poisson(const Grid& grid, const std::vector<double>& rhs,
                                std::vector<double>& x, double tolerance,
                                int max_iterations) {
    // Solve -Delta x = rhs in the zero-mean subspace. The periodic Laplacian's
    // constant null mode is removed explicitly.
    std::vector<double> b = rhs;
    remove_mean(b);
    std::fill(x.begin(), x.end(), 0.0);

    std::vector<double> r = b;
    std::vector<double> p = r;
    std::vector<double> ap(grid.cells, 0.0);
    const double b2 = dot(b, b);
    if (b2 <= std::numeric_limits<double>::min()) {
        return {};
    }

    double r2 = b2;
    CgResult result;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        apply_negative_laplacian(grid, p, ap);
        const double denominator = dot(p, ap);
        if (!(denominator > 0.0)) {
            result.iterations = iteration;
            result.relative_residual = std::sqrt(r2 / b2);
            result.converged = false;
            return result;
        }
        const double alpha = r2 / denominator;
        for (std::size_t q = 0; q < grid.cells; ++q) {
            x[q] += alpha * p[q];
            r[q] -= alpha * ap[q];
        }
        // Roundoff can slowly reintroduce the null mode.
        remove_mean(r);
        const double next_r2 = dot(r, r);
        result.iterations = iteration;
        result.relative_residual = std::sqrt(next_r2 / b2);
        if (result.relative_residual <= tolerance) {
            remove_mean(x);
            result.converged = true;
            return result;
        }
        const double beta = next_r2 / r2;
        for (std::size_t q = 0; q < grid.cells; ++q) {
            p[q] = r[q] + beta * p[q];
        }
        r2 = next_r2;
    }
    remove_mean(x);
    result.converged = false;
    return result;
}

void divergence_forward(const Grid& grid, const Field& velocity, std::vector<double>& result) {
    const double inv_h = 1.0 / grid.h;
    for (int k = 0; k < grid.n; ++k) {
        for (int j = 0; j < grid.n; ++j) {
            for (int i = 0; i < grid.n; ++i) {
                const auto q = grid.at(i, j, k);
                result[q] = ((velocity.c[0][grid.at(i + 1, j, k)] - velocity.c[0][q]) +
                             (velocity.c[1][grid.at(i, j + 1, k)] - velocity.c[1][q]) +
                             (velocity.c[2][grid.at(i, j, k + 1)] - velocity.c[2][q])) *
                            inv_h;
            }
        }
    }
}

CgResult project(const Grid& grid, Field& velocity, double scale, double tolerance,
                 int max_iterations) {
    std::vector<double> div(grid.cells, 0.0);
    divergence_forward(grid, velocity, div);
    std::vector<double> rhs(grid.cells, 0.0);
    for (std::size_t q = 0; q < grid.cells; ++q) {
        rhs[q] = -div[q] / scale;
    }
    std::vector<double> pressure(grid.cells, 0.0);
    const CgResult cg = solve_periodic_poisson(grid, rhs, pressure, tolerance, max_iterations);

    const double factor = scale / grid.h;
    for (int k = 0; k < grid.n; ++k) {
        for (int j = 0; j < grid.n; ++j) {
            for (int i = 0; i < grid.n; ++i) {
                const auto q = grid.at(i, j, k);
                velocity.c[0][q] -= factor * (pressure[q] - pressure[grid.at(i - 1, j, k)]);
                velocity.c[1][q] -= factor * (pressure[q] - pressure[grid.at(i, j - 1, k)]);
                velocity.c[2][q] -= factor * (pressure[q] - pressure[grid.at(i, j, k - 1)]);
            }
        }
    }
    return cg;
}

[[nodiscard]] double central_derivative(const Grid& grid, const std::vector<double>& f,
                                        int i, int j, int k, int direction) {
    if (direction == 0) {
        return (f[grid.at(i + 1, j, k)] - f[grid.at(i - 1, j, k)]) / (2.0 * grid.h);
    }
    if (direction == 1) {
        return (f[grid.at(i, j + 1, k)] - f[grid.at(i, j - 1, k)]) / (2.0 * grid.h);
    }
    return (f[grid.at(i, j, k + 1)] - f[grid.at(i, j, k - 1)]) / (2.0 * grid.h);
}

void compute_rhs(const Grid& grid, const Field& velocity, double viscosity, Field& rhs) {
    const double inv_h2 = 1.0 / (grid.h * grid.h);
    for (int k = 0; k < grid.n; ++k) {
        for (int j = 0; j < grid.n; ++j) {
            for (int i = 0; i < grid.n; ++i) {
                const auto q = grid.at(i, j, k);
                for (int component = 0; component < 3; ++component) {
                    const auto& f = velocity.c[static_cast<std::size_t>(component)];
                    double advection = 0.0;
                    for (int direction = 0; direction < 3; ++direction) {
                        const auto& transporting = velocity.c[static_cast<std::size_t>(direction)];
                        const double advective = transporting[q] *
                            central_derivative(grid, f, i, j, k, direction);

                        int ip = i;
                        int im = i;
                        int jp = j;
                        int jm = j;
                        int kp = k;
                        int km = k;
                        if (direction == 0) {
                            ++ip;
                            --im;
                        } else if (direction == 1) {
                            ++jp;
                            --jm;
                        } else {
                            ++kp;
                            --km;
                        }
                        const auto plus = grid.at(ip, jp, kp);
                        const auto minus = grid.at(im, jm, km);
                        const double conservative =
                            (transporting[plus] * f[plus] - transporting[minus] * f[minus]) /
                            (2.0 * grid.h);
                        advection += 0.5 * (advective + conservative);
                    }
                    const double laplacian =
                        (f[grid.at(i + 1, j, k)] + f[grid.at(i - 1, j, k)] +
                         f[grid.at(i, j + 1, k)] + f[grid.at(i, j - 1, k)] +
                         f[grid.at(i, j, k + 1)] + f[grid.at(i, j, k - 1)] - 6.0 * f[q]) *
                        inv_h2;
                    rhs.c[static_cast<std::size_t>(component)][q] =
                        -advection + viscosity * laplacian;
                }
            }
        }
    }
}

struct Diagnostics {
    double energy = 0.0;
    double enstrophy = 0.0;
    double gradient_squared = 0.0;
    double velocity_linf = 0.0;
    double vorticity_linf = 0.0;
    double divergence_l2 = 0.0;
};

Diagnostics diagnose(const Grid& grid, const Field& velocity) {
    Diagnostics d;
    double velocity2_sum = 0.0;
    double vorticity2_sum = 0.0;
    double gradient2_sum = 0.0;
    double divergence2_sum = 0.0;
    const double inv_h = 1.0 / grid.h;

    for (int k = 0; k < grid.n; ++k) {
        for (int j = 0; j < grid.n; ++j) {
            for (int i = 0; i < grid.n; ++i) {
                const auto q = grid.at(i, j, k);
                double speed2 = 0.0;
                for (int component = 0; component < 3; ++component) {
                    const auto& f = velocity.c[static_cast<std::size_t>(component)];
                    speed2 += f[q] * f[q];
                    const double dx = (f[grid.at(i + 1, j, k)] - f[q]) * inv_h;
                    const double dy = (f[grid.at(i, j + 1, k)] - f[q]) * inv_h;
                    const double dz = (f[grid.at(i, j, k + 1)] - f[q]) * inv_h;
                    gradient2_sum += dx * dx + dy * dy + dz * dz;
                }

                const double omega_x = central_derivative(grid, velocity.c[2], i, j, k, 1) -
                                       central_derivative(grid, velocity.c[1], i, j, k, 2);
                const double omega_y = central_derivative(grid, velocity.c[0], i, j, k, 2) -
                                       central_derivative(grid, velocity.c[2], i, j, k, 0);
                const double omega_z = central_derivative(grid, velocity.c[1], i, j, k, 0) -
                                       central_derivative(grid, velocity.c[0], i, j, k, 1);
                const double omega2 = omega_x * omega_x + omega_y * omega_y + omega_z * omega_z;
                const double div = ((velocity.c[0][grid.at(i + 1, j, k)] - velocity.c[0][q]) +
                                    (velocity.c[1][grid.at(i, j + 1, k)] - velocity.c[1][q]) +
                                    (velocity.c[2][grid.at(i, j, k + 1)] - velocity.c[2][q])) *
                                   inv_h;
                velocity2_sum += speed2;
                vorticity2_sum += omega2;
                divergence2_sum += div * div;
                d.velocity_linf = std::max(d.velocity_linf, std::sqrt(speed2));
                d.vorticity_linf = std::max(d.vorticity_linf, std::sqrt(omega2));
            }
        }
    }
    const double inverse_cells = 1.0 / static_cast<double>(grid.cells);
    d.energy = 0.5 * velocity2_sum * inverse_cells;
    d.enstrophy = 0.5 * vorticity2_sum * inverse_cells;
    d.gradient_squared = gradient2_sum * inverse_cells;
    d.divergence_l2 = std::sqrt(divergence2_sum * inverse_cells);
    return d;
}

[[nodiscard]] double field_mean_square(const Grid& grid, const Field& field) {
    double sum = 0.0;
    for (const auto& component : field.c) {
        sum += dot(component, component);
    }
    return sum / static_cast<double>(grid.cells);
}

struct Config {
    int n = 16;
    double viscosity = 0.02;
    double final_time = 0.25;
    double cfl = 0.35;
    double max_dt = 0.005;
    double amplitude = 1.0;
    std::string initial_condition = "taylor-green";
    std::uint64_t seed = 1;
    int modes = 8;
    int output_every = 10;
    double poisson_tolerance = 1e-10;
    int poisson_iterations = 1000;
    bool quiet = false;
    std::string csv_path;
};

void initialize_taylor_green(const Grid& grid, Field& velocity, double amplitude) {
    for (int k = 0; k < grid.n; ++k) {
        for (int j = 0; j < grid.n; ++j) {
            for (int i = 0; i < grid.n; ++i) {
                const double x = (static_cast<double>(i) + 0.5) * grid.h;
                const double y = (static_cast<double>(j) + 0.5) * grid.h;
                const double z = (static_cast<double>(k) + 0.5) * grid.h;
                const auto q = grid.at(i, j, k);
                velocity.c[0][q] = amplitude * std::sin(x) * std::cos(y) * std::cos(z);
                velocity.c[1][q] = -amplitude * std::cos(x) * std::sin(y) * std::cos(z);
                velocity.c[2][q] = 0.0;
            }
        }
    }
}

void initialize_random_modes(const Grid& grid, Field& velocity, double amplitude,
                             std::uint64_t seed, int mode_count) {
    if (mode_count < 1) {
        throw std::invalid_argument("--modes must be positive");
    }
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<int> wave_number(-3, 3);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> phase_distribution(0.0, 2.0 * pi);
    const double mode_scale = amplitude / std::sqrt(static_cast<double>(mode_count));

    for (int mode = 0; mode < mode_count; ++mode) {
        std::array<int, 3> wave{};
        do {
            for (int direction = 0; direction < 3; ++direction) {
                wave[static_cast<std::size_t>(direction)] = wave_number(generator);
            }
        } while (wave[0] == 0 && wave[1] == 0 && wave[2] == 0);

        std::array<double, 3> polarization{normal(generator), normal(generator), normal(generator)};
        const double wave2 = static_cast<double>(wave[0] * wave[0] + wave[1] * wave[1] +
                                                 wave[2] * wave[2]);
        const double projection =
            (polarization[0] * static_cast<double>(wave[0]) +
             polarization[1] * static_cast<double>(wave[1]) +
             polarization[2] * static_cast<double>(wave[2])) /
            wave2;
        double polarization2 = 0.0;
        for (int direction = 0; direction < 3; ++direction) {
            polarization[static_cast<std::size_t>(direction)] -=
                projection * static_cast<double>(wave[static_cast<std::size_t>(direction)]);
            polarization2 += polarization[static_cast<std::size_t>(direction)] *
                             polarization[static_cast<std::size_t>(direction)];
        }
        if (polarization2 < 1e-14) {
            --mode;
            continue;
        }
        const double normalization = mode_scale / std::sqrt(polarization2);
        const double phase = phase_distribution(generator);
        for (int k = 0; k < grid.n; ++k) {
            for (int j = 0; j < grid.n; ++j) {
                for (int i = 0; i < grid.n; ++i) {
                    const double x = (static_cast<double>(i) + 0.5) * grid.h;
                    const double y = (static_cast<double>(j) + 0.5) * grid.h;
                    const double z = (static_cast<double>(k) + 0.5) * grid.h;
                    const double angle = static_cast<double>(wave[0]) * x +
                                         static_cast<double>(wave[1]) * y +
                                         static_cast<double>(wave[2]) * z + phase;
                    const auto q = grid.at(i, j, k);
                    for (int component = 0; component < 3; ++component) {
                        velocity.c[static_cast<std::size_t>(component)][q] +=
                            normalization * polarization[static_cast<std::size_t>(component)] *
                            std::cos(angle);
                    }
                }
            }
        }
    }
}

struct RunSummary {
    Config config;
    int steps = 0;
    double simulated_time = 0.0;
    Diagnostics initial;
    Diagnostics final;
    double peak_energy = 0.0;
    double peak_enstrophy = 0.0;
    double peak_velocity_linf = 0.0;
    double peak_vorticity_linf = 0.0;
    double prodi_serrin_integral = 0.0;
    double bkm_integral = 0.0;
    double max_discrete_energy_constant = 0.0;
    double max_energy_lemma_excess = 0.0;
    double max_enstrophy_growth_ratio = 0.0;
    int max_poisson_iterations = 0;
    bool poisson_converged = true;
    bool finite = true;
};

void print_diagnostic_header(std::ostream& out) {
    out << "step,time,dt,energy,enstrophy,velocity_linf,vorticity_linf,divergence_l2,"
           "energy_lemma_C,poisson_iterations\n";
}

void print_diagnostic_row(std::ostream& out, int step, double time, double dt,
                          const Diagnostics& d, double energy_constant, int poisson_iterations) {
    out << step << ',' << std::setprecision(12) << time << ',' << dt << ',' << d.energy << ','
        << d.enstrophy << ',' << d.velocity_linf << ',' << d.vorticity_linf << ','
        << d.divergence_l2 << ',' << energy_constant << ',' << poisson_iterations << '\n';
}

RunSummary simulate(const Config& config) {
    const Grid grid(config.n);
    Field velocity(grid.cells);
    if (config.initial_condition == "taylor-green") {
        initialize_taylor_green(grid, velocity, config.amplitude);
    } else if (config.initial_condition == "random") {
        initialize_random_modes(grid, velocity, config.amplitude, config.seed, config.modes);
    } else {
        throw std::invalid_argument("unknown initial condition: " + config.initial_condition);
    }

    const CgResult initial_projection = project(grid, velocity, 1.0, config.poisson_tolerance,
                                                config.poisson_iterations);
    RunSummary summary;
    summary.config = config;
    summary.poisson_converged = initial_projection.converged;
    summary.max_poisson_iterations = initial_projection.iterations;
    summary.initial = diagnose(grid, velocity);
    summary.final = summary.initial;
    summary.peak_energy = summary.initial.energy;
    summary.peak_enstrophy = summary.initial.enstrophy;
    summary.peak_velocity_linf = summary.initial.velocity_linf;
    summary.peak_vorticity_linf = summary.initial.vorticity_linf;

    std::ofstream csv;
    if (!config.csv_path.empty()) {
        csv.open(config.csv_path);
        if (!csv) {
            throw std::runtime_error("cannot open CSV output: " + config.csv_path);
        }
        print_diagnostic_header(csv);
        print_diagnostic_row(csv, 0, 0.0, 0.0, summary.initial, 0.0,
                             initial_projection.iterations);
    }
    if (!config.quiet) {
        print_diagnostic_header(std::cout);
        print_diagnostic_row(std::cout, 0, 0.0, 0.0, summary.initial, 0.0,
                             initial_projection.iterations);
    }

    Field rhs(grid.cells);
    double time = 0.0;
    int step = 0;
    while (time < config.final_time) {
        const Diagnostics before = diagnose(grid, velocity);
        double dt_advection = config.max_dt;
        if (before.velocity_linf > 1e-14) {
            dt_advection = config.cfl * grid.h / before.velocity_linf;
        }
        double dt_diffusion = config.max_dt;
        if (config.viscosity > 0.0) {
            dt_diffusion = 0.9 * grid.h * grid.h / (6.0 * config.viscosity);
        }
        const double dt = std::min({config.max_dt, dt_advection, dt_diffusion,
                                    config.final_time - time});
        if (!(dt > 0.0) || !std::isfinite(dt)) {
            summary.finite = false;
            break;
        }

        compute_rhs(grid, velocity, config.viscosity, rhs);
        const double rhs_mean_square = field_mean_square(grid, rhs);
        for (int component = 0; component < 3; ++component) {
            for (std::size_t q = 0; q < grid.cells; ++q) {
                velocity.c[static_cast<std::size_t>(component)][q] +=
                    dt * rhs.c[static_cast<std::size_t>(component)][q];
            }
        }
        const CgResult cg = project(grid, velocity, dt, config.poisson_tolerance,
                                    config.poisson_iterations);
        const Diagnostics after = diagnose(grid, velocity);

        // Discrete one-step lemma:
        // E_{n+1} - E_n + nu*dt*||grad u_n||_2^2
        //     <= 0.5*dt^2*||-N(u_n)+nu*Delta u_n||_2^2.
        // The skew advection has zero discrete energy pairing and projection is
        // an L2 contraction. Hence C <= 1 up to Poisson/roundoff error.
        const double lemma_left = after.energy - before.energy +
                                  config.viscosity * dt * before.gradient_squared;
        const double lemma_right = 0.5 * dt * dt * rhs_mean_square;
        double energy_constant = 0.0;
        if (lemma_right > 1e-30) {
            energy_constant = lemma_left / lemma_right;
        }
        summary.max_discrete_energy_constant =
            std::max(summary.max_discrete_energy_constant, energy_constant);
        summary.max_energy_lemma_excess =
            std::max(summary.max_energy_lemma_excess, lemma_left - lemma_right);

        const double enstrophy_growth = (after.enstrophy - before.enstrophy) / dt;
        if (enstrophy_growth > 0.0 && before.enstrophy > 1e-14 && config.viscosity > 0.0) {
            const double ratio = enstrophy_growth * std::pow(config.viscosity, 3.0) /
                                 std::pow(before.enstrophy, 3.0);
            summary.max_enstrophy_growth_ratio =
                std::max(summary.max_enstrophy_growth_ratio, ratio);
        }

        summary.prodi_serrin_integral +=
            0.5 * dt * (before.velocity_linf * before.velocity_linf +
                        after.velocity_linf * after.velocity_linf);
        summary.bkm_integral +=
            0.5 * dt * (before.vorticity_linf + after.vorticity_linf);
        summary.poisson_converged = summary.poisson_converged && cg.converged;
        summary.max_poisson_iterations = std::max(summary.max_poisson_iterations, cg.iterations);
        summary.peak_energy = std::max(summary.peak_energy, after.energy);
        summary.peak_enstrophy = std::max(summary.peak_enstrophy, after.enstrophy);
        summary.peak_velocity_linf = std::max(summary.peak_velocity_linf, after.velocity_linf);
        summary.peak_vorticity_linf =
            std::max(summary.peak_vorticity_linf, after.vorticity_linf);

        ++step;
        time += dt;
        summary.steps = step;
        summary.simulated_time = time;
        summary.final = after;
        summary.finite = summary.finite && std::isfinite(after.energy) &&
                         std::isfinite(after.enstrophy) && std::isfinite(after.divergence_l2);

        const bool output_now = step % config.output_every == 0 || time >= config.final_time ||
                                !summary.finite;
        if (output_now && !config.quiet) {
            print_diagnostic_row(std::cout, step, time, dt, after, energy_constant, cg.iterations);
        }
        if (csv && output_now) {
            print_diagnostic_row(csv, step, time, dt, after, energy_constant, cg.iterations);
        }
        if (!summary.finite || !cg.converged) {
            break;
        }
    }
    return summary;
}

void print_summary(const RunSummary& summary, std::ostream& out) {
    out << std::setprecision(10)
        << "\nRun summary\n"
        << "  grid:                       " << summary.config.n << "^3\n"
        << "  steps / time:               " << summary.steps << " / "
        << summary.simulated_time << "\n"
        << "  peak energy:                " << summary.peak_energy << "\n"
        << "  peak enstrophy:             " << summary.peak_enstrophy << "\n"
        << "  peak ||u||_inf:             " << summary.peak_velocity_linf << "\n"
        << "  peak ||omega||_inf:         " << summary.peak_vorticity_linf << "\n"
        << "  integral ||u||_inf^2 dt:    " << summary.prodi_serrin_integral << "\n"
        << "  integral ||omega||_inf dt:  " << summary.bkm_integral << "\n"
        << "  max discrete lemma C:       " << summary.max_discrete_energy_constant << " (target <= 1)\n"
        << "  max lemma absolute excess:  " << summary.max_energy_lemma_excess << "\n"
        << "  max normalized dZ/dt ratio: " << summary.max_enstrophy_growth_ratio << "\n"
        << "  final divergence L2:        " << summary.final.divergence_l2 << "\n"
        << "  Poisson solver:             "
        << (summary.poisson_converged ? "converged" : "FAILED") << "\n"
        << "  finite values:              " << (summary.finite ? "yes" : "NO") << "\n";
}

[[nodiscard]] std::vector<int> parse_int_list(const std::string& text) {
    std::vector<int> values;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        values.push_back(std::stoi(item));
    }
    if (values.empty()) {
        throw std::invalid_argument("empty integer list");
    }
    return values;
}

[[nodiscard]] std::vector<std::uint64_t> parse_seed_list(const std::string& text) {
    std::vector<std::uint64_t> values;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        values.push_back(static_cast<std::uint64_t>(std::stoull(item)));
    }
    if (values.empty()) {
        throw std::invalid_argument("empty seed list");
    }
    return values;
}

struct ParsedOptions {
    Config config;
    std::vector<int> grids{8, 12, 16};
    std::vector<std::uint64_t> seeds{1, 2, 3};
    std::string report_path;
};

[[nodiscard]] ParsedOptions parse_options(int argc, char** argv, int first) {
    ParsedOptions result;
    auto value_after = [&](int& index, std::string_view option) -> std::string {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value after " + std::string(option));
        }
        ++index;
        return argv[index];
    };
    for (int i = first; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--n") {
            result.config.n = std::stoi(value_after(i, option));
        } else if (option == "--nu") {
            result.config.viscosity = std::stod(value_after(i, option));
        } else if (option == "--time") {
            result.config.final_time = std::stod(value_after(i, option));
        } else if (option == "--cfl") {
            result.config.cfl = std::stod(value_after(i, option));
        } else if (option == "--max-dt") {
            result.config.max_dt = std::stod(value_after(i, option));
        } else if (option == "--amplitude") {
            result.config.amplitude = std::stod(value_after(i, option));
        } else if (option == "--init") {
            result.config.initial_condition = value_after(i, option);
        } else if (option == "--seed") {
            result.config.seed = static_cast<std::uint64_t>(std::stoull(value_after(i, option)));
        } else if (option == "--modes") {
            result.config.modes = std::stoi(value_after(i, option));
        } else if (option == "--output-every") {
            result.config.output_every = std::stoi(value_after(i, option));
        } else if (option == "--poisson-tol") {
            result.config.poisson_tolerance = std::stod(value_after(i, option));
        } else if (option == "--poisson-iters") {
            result.config.poisson_iterations = std::stoi(value_after(i, option));
        } else if (option == "--csv") {
            result.config.csv_path = value_after(i, option);
        } else if (option == "--quiet") {
            result.config.quiet = true;
        } else if (option == "--grids") {
            result.grids = parse_int_list(value_after(i, option));
        } else if (option == "--seeds") {
            result.seeds = parse_seed_list(value_after(i, option));
        } else if (option == "--report") {
            result.report_path = value_after(i, option);
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    if (result.config.viscosity < 0.0 || result.config.final_time <= 0.0 ||
        result.config.cfl <= 0.0 || result.config.max_dt <= 0.0 ||
        result.config.output_every <= 0 || result.config.poisson_tolerance <= 0.0 ||
        result.config.poisson_iterations <= 0) {
        throw std::invalid_argument("numeric parameters are outside their valid range");
    }
    return result;
}

[[nodiscard]] double log_slope(int n0, double y0, int n1, double y1) {
    if (y0 <= 0.0 || y1 <= 0.0 || n0 == n1) {
        return 0.0;
    }
    return std::log(y1 / y0) / std::log(static_cast<double>(n1) / static_cast<double>(n0));
}

void write_mining_report(std::ostream& out, const std::vector<RunSummary>& runs,
                         const std::vector<int>& grids, const std::vector<std::uint64_t>& seeds) {
    out << "# Navier-Stokes lemma-mining report\n\n"
        << "This is a numerical falsification report, not a proof of 3D regularity.\n\n"
        << "| N | seed | peak Z | peak omega_inf | int u_inf^2 | energy C | div L2 | status |\n"
        << "|---:|---:|---:|---:|---:|---:|---:|:---|\n";
    for (const auto& run : runs) {
        out << '|' << run.config.n << '|' << run.config.seed << '|' << std::setprecision(8)
            << run.peak_enstrophy << '|' << run.peak_vorticity_linf << '|'
            << run.prodi_serrin_integral << '|' << run.max_discrete_energy_constant << '|'
            << run.final.divergence_l2 << '|'
            << (run.finite && run.poisson_converged && run.max_discrete_energy_constant <= 1.00001
                    ? "pass"
                    : "inspect")
            << "|\n";
    }

    out << "\n## Resolution trend\n\n"
        << "The useful continuous lemma candidate would be a bound uniform in N for one of the "
           "conditional regularity quantities. A positive fitted power signals that the current "
           "data do not yet support such a bound.\n\n"
        << "| seed | slope peak Z | slope peak omega_inf | slope int u_inf^2 | verdict |\n"
        << "|---:|---:|---:|---:|:---|\n";

    for (const auto seed : seeds) {
        std::vector<const RunSummary*> selected;
        for (const auto& run : runs) {
            if (run.config.seed == seed) {
                selected.push_back(&run);
            }
        }
        std::sort(selected.begin(), selected.end(), [](const auto* a, const auto* b) {
            return a->config.n < b->config.n;
        });
        if (selected.size() < 2) {
            continue;
        }
        const auto& low = *selected.front();
        const auto& high = *selected.back();
        const double slope_z = log_slope(low.config.n, low.peak_enstrophy,
                                         high.config.n, high.peak_enstrophy);
        const double slope_omega = log_slope(low.config.n, low.peak_vorticity_linf,
                                             high.config.n, high.peak_vorticity_linf);
        const double slope_ps = log_slope(low.config.n, low.prodi_serrin_integral,
                                          high.config.n, high.prodi_serrin_integral);
        const bool bounded_candidate = slope_z < 0.25 && slope_omega < 0.25 && slope_ps < 0.25;
        out << '|' << seed << '|' << slope_z << '|' << slope_omega << '|' << slope_ps << '|'
            << (bounded_candidate ? "candidate survives" : "needs refinement / rejected") << "|\n";
    }

    out << "\n## What has and has not been established\n\n"
        << "The tested discrete lemma is\n\n"
        << "```text\n"
        << "E[n+1] - E[n] + nu*dt*||grad u[n]||_2^2\n"
        << "    <= (dt^2/2)*||-N_h(u[n]) + nu*Delta_h u[n]||_2^2.\n"
        << "```\n\n"
        << "For this finite-dimensional scheme it follows from skew-symmetry of `N_h` and the "
           "L2-contractivity of the pressure projection. The experiment checks the implementation. "
           "Passing it does not give an N-independent bound on enstrophy or a strong-solution norm. "
           "That missing uniform bound is the analytical gap to the Millennium problem.\n\n"
        << "Grids tested: ";
    for (std::size_t q = 0; q < grids.size(); ++q) {
        out << (q == 0 ? "" : ", ") << grids[q] << "^3";
    }
    out << ".\n";
}

int run_miner(const ParsedOptions& options) {
    std::vector<RunSummary> runs;
    std::cout << "N,seed,peak_enstrophy,peak_omega_inf,prodi_serrin,energy_C,div_l2,status\n";
    for (const int n : options.grids) {
        for (const auto seed : options.seeds) {
            Config config = options.config;
            config.n = n;
            config.seed = seed;
            config.initial_condition = "random";
            config.quiet = true;
            config.csv_path.clear();
            RunSummary run = simulate(config);
            std::cout << n << ',' << seed << ',' << std::setprecision(10) << run.peak_enstrophy
                      << ',' << run.peak_vorticity_linf << ',' << run.prodi_serrin_integral << ','
                      << run.max_discrete_energy_constant << ',' << run.final.divergence_l2 << ','
                      << (run.finite && run.poisson_converged ? "ok" : "failed") << '\n';
            runs.push_back(std::move(run));
        }
    }

    if (!options.report_path.empty()) {
        std::ofstream report(options.report_path);
        if (!report) {
            throw std::runtime_error("cannot open report output: " + options.report_path);
        }
        write_mining_report(report, runs, options.grids, options.seeds);
        std::cout << "report written to " << options.report_path << '\n';
    } else {
        std::cout << '\n';
        write_mining_report(std::cout, runs, options.grids, options.seeds);
    }
    const bool all_good = std::all_of(runs.begin(), runs.end(), [](const RunSummary& run) {
        return run.finite && run.poisson_converged &&
               run.max_discrete_energy_constant <= 1.00001 && run.final.divergence_l2 < 1e-7;
    });
    return all_good ? 0 : 2;
}

int self_test() {
    bool okay = true;
    {
        const Grid grid(8);
        Field field(grid.cells);
        initialize_random_modes(grid, field, 1.0, 42, 5);
        const Diagnostics before = diagnose(grid, field);
        const CgResult cg = project(grid, field, 1.0, 1e-12, 1000);
        const Diagnostics after = diagnose(grid, field);
        const bool pass = cg.converged && after.divergence_l2 < 1e-10 &&
                          after.divergence_l2 < before.divergence_l2 * 1e-8;
        std::cout << "projection test: " << (pass ? "PASS" : "FAIL")
                  << " (div " << before.divergence_l2 << " -> " << after.divergence_l2 << ")\n";
        okay = okay && pass;
    }
    {
        Config config;
        config.n = 8;
        config.final_time = 0.03;
        config.max_dt = 0.003;
        config.quiet = true;
        config.poisson_tolerance = 1e-12;
        const RunSummary run = simulate(config);
        const bool pass = run.finite && run.poisson_converged &&
                          run.max_discrete_energy_constant <= 1.00001 &&
                          run.final.divergence_l2 < 1e-9 &&
                          run.final.energy <= run.initial.energy + 1e-3;
        std::cout << "time-step/lemma test: " << (pass ? "PASS" : "FAIL")
                  << " (C=" << run.max_discrete_energy_constant
                  << ", div=" << run.final.divergence_l2 << ")\n";
        okay = okay && pass;
    }
    okay = lemma::self_test(std::cout) && okay;
    return okay ? 0 : 1;
}

void print_help(std::ostream& out) {
    out << "Navier-Stokes lemma laboratory (periodic 3D finite differences)\n\n"
        << "Usage:\n"
        << "  navier_stokes_lab simulate [options]\n"
        << "  navier_stokes_lab mine [options]\n"
        << "  navier_stokes_lab lemma [options]\n"
        << "  navier_stokes_lab adversary [options]\n"
        << "  navier_stokes_lab family [options]\n"
        << "  navier_stokes_lab helical-adversary [options]\n"
        << "  navier_stokes_lab helical-cutoff-scan [options]\n"
        << "  navier_stokes_lab orthogonal-triad-certificate [options]\n"
        << "  navier_stokes_lab local-signature-certificate [options]\n"
        << "  navier_stokes_lab local-signature-adversary [options]\n"
        << "  navier_stokes_lab local-signature-factor [options]\n"
        << "  navier_stokes_lab local-signature-gradient [options]\n"
        << "  navier_stokes_lab local-signature-trajectory [options]\n"
        << "  navier_stokes_lab local-closure-adversary [options]\n"
        << "  navier_stokes_lab local-sld-ansatz [options]\n"
        << "  navier_stokes_lab local-sld-trajectory-ansatz [options]\n"
        << "  navier_stokes_lab local-sld-krylov-ansatz [options]\n"
        << "  navier_stokes_lab local-sld-response-hierarchy [options]\n"
        << "  navier_stokes_lab local-sld-block [options]\n"
        << "  navier_stokes_lab shifted-density [options]\n"
        << "  navier_stokes_lab self-test\n\n"
        << "Simulation options:\n"
        << "  --n N                 grid N^3 (default 16)\n"
        << "  --nu VALUE            viscosity (default 0.02)\n"
        << "  --time VALUE          final time (default 0.25)\n"
        << "  --cfl VALUE           advective CFL (default 0.35)\n"
        << "  --max-dt VALUE        maximum time step (default 0.005)\n"
        << "  --init NAME           taylor-green or random\n"
        << "  --amplitude VALUE     initial amplitude\n"
        << "  --seed N              random seed\n"
        << "  --modes N             number of random Fourier modes\n"
        << "  --output-every N      diagnostic interval\n"
        << "  --csv PATH            write diagnostics to CSV\n"
        << "  --quiet               suppress per-step table\n\n"
        << "Mining options:\n"
        << "  --grids A,B,C         resolution ladder (default 8,12,16)\n"
        << "  --seeds A,B,C         ensemble seeds (default 1,2,3)\n"
        << "  --report PATH         write Markdown report\n";
    out << '\n';
    lemma::LemmaCli::print_help(out);
    out << '\n';
    lemma::LemmaCli::print_adversary_help(out);
    out << '\n';
    lemma::LemmaCli::print_family_help(out);
    out << '\n';
    lemma::HelicalAdversaryCli::print_help(out);
    out << '\n';
    lemma::HelicalCutoffScan::print_help(out);
    out << '\n';
    lemma::OrthogonalTriadCli::print_help(out);
    out << '\n';
    lemma::LocalSignatureCli::print_help(out);
    out << '\n';
    lemma::LocalSignatureAdversaryCli::print_help(out);
    out << '\n';
    lemma::LocalSignatureFactorAdversaryCli::print_help(out);
    out << '\n';
    lemma::LocalSignatureGradientCli::print_help(out);
    out << '\n';
    lemma::LocalSignatureTrajectoryCli::print_help(out);
    out << '\n';
    lemma::LocalQuarticClosureCli::print_help(out);
    out << '\n';
    lemma::LocalSldCyclicAnsatzCli::print_help(out);
    out << '\n';
    lemma::LocalSldCyclicTrajectoryCli::print_help(out);
    out << '\n';
    lemma::LocalSldCyclicKrylovCli::print_help(out);
    out << '\n';
    lemma::LocalSldResponseHierarchyCli::print_help(out);
    out << '\n';
    lemma::LocalSldSignatureBlockCli::print_help(out);
    out << '\n';
    lemma::ShiftedCriticalDensityCli::print_help(out);
    out << '\n';
    lemma::StateAnalysisCli::print_help(out);
    out << '\n';
    lemma::StateFamilyAnalysisCli::print_help(out);
}

}  // namespace ns

int main(int argc, char** argv) {
    try {
        if (argc < 2 || std::string_view(argv[1]) == "--help" ||
            std::string_view(argv[1]) == "-h") {
            ns::print_help(std::cout);
            return argc < 2 ? 1 : 0;
        }
        const std::string command = argv[1];
        if (command == "self-test") {
            return ns::self_test();
        }
        if (command == "lemma") {
            return lemma::run(
                lemma::LemmaCli::parse_options(argc, argv, 2), std::cout);
        }
        if (command == "adversary") {
            return lemma::run_adversary(
                lemma::LemmaCli::parse_adversary_options(argc, argv, 2), std::cout);
        }
        if (command == "family") {
            return lemma::run_family(
                lemma::LemmaCli::parse_family_options(argc, argv, 2), std::cout);
        }
        if (command == "helical-adversary") {
            return lemma::run_helical_adversary(
                lemma::HelicalAdversaryCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "helical-cutoff-scan") {
            return lemma::HelicalCutoffScan::run(
                lemma::HelicalCutoffScan::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "orthogonal-triad-certificate") {
            return lemma::OrthogonalTriadCli::run(
                lemma::OrthogonalTriadCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-signature-certificate") {
            return lemma::LocalSignatureCli::run(
                lemma::LocalSignatureCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-signature-adversary") {
            return lemma::LocalSignatureAdversaryCli::run(
                lemma::LocalSignatureAdversaryCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-signature-factor") {
            return lemma::LocalSignatureFactorAdversaryCli::run(
                lemma::LocalSignatureFactorAdversaryCli::parse(
                    argc, argv, 2),
                std::cout);
        }
        if (command == "local-signature-gradient") {
            return lemma::LocalSignatureGradientCli::run(
                lemma::LocalSignatureGradientCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-signature-trajectory") {
            return lemma::LocalSignatureTrajectoryCli::run(
                lemma::LocalSignatureTrajectoryCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-closure-adversary") {
            return lemma::LocalQuarticClosureCli::run(
                lemma::LocalQuarticClosureCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-sld-ansatz") {
            return lemma::LocalSldCyclicAnsatzCli::run(
                lemma::LocalSldCyclicAnsatzCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "local-sld-trajectory-ansatz") {
            return lemma::LocalSldCyclicTrajectoryCli::run(
                lemma::LocalSldCyclicTrajectoryCli::parse(
                    argc, argv, 2),
                std::cout);
        }
        if (command == "local-sld-krylov-ansatz") {
            return lemma::LocalSldCyclicKrylovCli::run(
                lemma::LocalSldCyclicKrylovCli::parse(
                    argc, argv, 2),
                std::cout);
        }
        if (command == "local-sld-response-hierarchy") {
            return lemma::LocalSldResponseHierarchyCli::run(
                lemma::LocalSldResponseHierarchyCli::parse(
                    argc, argv, 2),
                std::cout);
        }
        if (command == "local-sld-block") {
            return lemma::LocalSldSignatureBlockCli::run(
                lemma::LocalSldSignatureBlockCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "shifted-density") {
            return lemma::ShiftedCriticalDensityCli::run(
                lemma::ShiftedCriticalDensityCli::parse(argc, argv, 2),
                std::cout);
        }
        if (command == "state-analysis") {
            return lemma::run_state_analysis(
                lemma::StateAnalysisCli::parse(argc, argv, 2), std::cout);
        }
        if (command == "state-family-analysis") {
            return lemma::run_state_family_analysis(
                lemma::StateFamilyAnalysisCli::parse(argc, argv, 2),
                std::cout);
        }
        const ns::ParsedOptions options = ns::parse_options(argc, argv, 2);
        if (command == "simulate") {
            const ns::RunSummary summary = ns::simulate(options.config);
            ns::print_summary(summary, std::cout);
            return summary.finite && summary.poisson_converged ? 0 : 2;
        }
        if (command == "mine") {
            return ns::run_miner(options);
        }
        throw std::invalid_argument("unknown command: " + command);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
