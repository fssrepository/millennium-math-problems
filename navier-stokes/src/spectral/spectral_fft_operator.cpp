#include "spectral_fft_operator.hpp"

#include "fft3d.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using GridField = std::array<std::vector<SpectralComplex>, 3>;

struct PhysicalFields {
    GridField velocity;
    GridField vorticity;
    std::size_t grid_n = 0;
    std::size_t grid_cells = 0;
};

std::size_t wrapped(int wave_number, std::size_t grid_n) {
    const int n = static_cast<int>(grid_n);
    const int result = wave_number % n;
    return static_cast<std::size_t>(result < 0 ? result + n : result);
}

std::size_t grid_index(std::size_t x, std::size_t y, std::size_t z,
                       std::size_t grid_n) {
    return (z * grid_n + y) * grid_n + x;
}

std::size_t mode_cell(WaveVector wave, std::size_t grid_n) {
    return grid_index(wrapped(wave.x, grid_n), wrapped(wave.y, grid_n),
                      wrapped(wave.z, grid_n), grid_n);
}

ComplexVector curl_coefficient(WaveVector wave,
                               const ComplexVector& velocity) {
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    return {
        imaginary_unit *
            (static_cast<SpectralReal>(wave.y) * velocity[2] -
             static_cast<SpectralReal>(wave.z) * velocity[1]),
        imaginary_unit *
            (static_cast<SpectralReal>(wave.z) * velocity[0] -
             static_cast<SpectralReal>(wave.x) * velocity[2]),
        imaginary_unit *
            (static_cast<SpectralReal>(wave.x) * velocity[1] -
             static_cast<SpectralReal>(wave.y) * velocity[0])};
}

void transform_fields(GridField& first, GridField& second,
                      std::size_t grid_n, bool inverse,
                      int compute_threads) {
    const int transform_threads = std::min(6, compute_threads);
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for num_threads(transform_threads) schedule(static)
#endif
    for (int transform = 0; transform < 6; ++transform) {
        if (transform < 3) {
            fft3d::transform_3d(
                first[static_cast<std::size_t>(transform)], grid_n, inverse);
        } else {
            fft3d::transform_3d(
                second[static_cast<std::size_t>(transform - 3)], grid_n,
                inverse);
        }
    }
}

PhysicalFields to_physical(const SpectralState& state,
                           const SpectralIncrement& coefficients,
                           int compute_threads) {
    if (coefficients.size() != state.waves.size()) {
        throw std::invalid_argument("FFT coefficient layout mismatch");
    }
    PhysicalFields result;
    const int cutoff = SpectralStateOps::cutoff(state);
    result.grid_n = fft3d::next_power_of_two(
        static_cast<std::size_t>(3 * cutoff + 1));
    result.grid_cells =
        result.grid_n * result.grid_n * result.grid_n;
    for (GridField* field : {&result.velocity, &result.vorticity}) {
        for (auto& component : *field) {
            component.assign(result.grid_cells, {});
        }
    }
    const SpectralReal embedding_scale =
        static_cast<SpectralReal>(result.grid_cells);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const WaveVector wave = state.waves[mode];
        const std::size_t cell = mode_cell(wave, result.grid_n);
        const ComplexVector omega =
            curl_coefficient(wave, coefficients[mode]);
        for (std::size_t component = 0; component < 3; ++component) {
            result.velocity[component][cell] =
                embedding_scale * coefficients[mode][component];
            result.vorticity[component][cell] =
                embedding_scale * omega[component];
        }
    }
    transform_fields(result.velocity, result.vorticity, result.grid_n, true,
                     compute_threads);
    return result;
}

ComplexVector cross(const ComplexVector& left, const ComplexVector& right) {
    return {left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0]};
}

ComplexVector grid_value(const GridField& field, std::size_t cell) {
    return {field[0][cell], field[1][cell], field[2][cell]};
}

void set_grid_value(GridField& field, std::size_t cell,
                    const ComplexVector& value) {
    for (std::size_t component = 0; component < 3; ++component) {
        field[component][cell] = value[component];
    }
}

void transform_single_field(GridField& field, std::size_t grid_n,
                            bool inverse, int compute_threads) {
    const int transform_threads = std::min(3, compute_threads);
#ifdef NS_HAVE_OPENMP
#pragma omp parallel for num_threads(transform_threads) schedule(static)
#endif
    for (int component = 0; component < 3; ++component) {
        fft3d::transform_3d(field[static_cast<std::size_t>(component)],
                            grid_n, inverse);
    }
}

SpectralIncrement extract_projected(const SpectralState& state,
                                    const GridField& fourier_grid,
                                    std::size_t grid_n,
                                    SpectralReal scale) {
    SpectralIncrement result(state.waves.size());
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const WaveVector wave = state.waves[mode];
        ComplexVector coefficient =
            grid_value(fourier_grid, mode_cell(wave, grid_n));
        for (SpectralComplex& component : coefficient) {
            component *= scale;
        }
        result[mode] = project_divergence_free(wave, coefficient);
    }
    return result;
}

GridField cross_field(const PhysicalFields& left,
                      const PhysicalFields* tangent) {
    GridField result;
    for (auto& component : result) {
        component.resize(left.grid_cells);
    }
    for (std::size_t cell = 0; cell < left.grid_cells; ++cell) {
        ComplexVector value =
            cross(grid_value(left.vorticity, cell),
                  grid_value(left.velocity, cell));
        if (tangent != nullptr) {
            const ComplexVector first =
                cross(grid_value(tangent->vorticity, cell),
                      grid_value(left.velocity, cell));
            const ComplexVector second =
                cross(grid_value(left.vorticity, cell),
                      grid_value(tangent->velocity, cell));
            for (std::size_t component = 0; component < 3; ++component) {
                value[component] = first[component] + second[component];
            }
        }
        set_grid_value(result, cell, value);
    }
    return result;
}

}  // namespace

SpectralIncrement SpectralFftOperator::advection(
    const SpectralState& state, int compute_threads) {
    const PhysicalFields fields =
        to_physical(state, state.velocity, compute_threads);
    GridField nonlinear = cross_field(fields, nullptr);
    transform_single_field(nonlinear, fields.grid_n, false, compute_threads);
    return extract_projected(
        state, nonlinear, fields.grid_n,
        1.0L / static_cast<SpectralReal>(fields.grid_cells));
}

SpectralIncrement SpectralFftOperator::advection_jvp(
    const SpectralState& state, const SpectralIncrement& direction,
    int compute_threads) {
    const PhysicalFields fields =
        to_physical(state, state.velocity, compute_threads);
    const PhysicalFields tangent =
        to_physical(state, direction, compute_threads);
    GridField nonlinear_tangent = cross_field(fields, &tangent);
    transform_single_field(nonlinear_tangent, fields.grid_n, false,
                           compute_threads);
    return extract_projected(
        state, nonlinear_tangent, fields.grid_n,
        1.0L / static_cast<SpectralReal>(fields.grid_cells));
}

SpectralIncrement SpectralFftOperator::advection_vjp(
    const SpectralState& state,
    const SpectralIncrement& output_cotangent, int compute_threads) {
    if (output_cotangent.size() != state.waves.size()) {
        throw std::invalid_argument("FFT cotangent layout mismatch");
    }
    const PhysicalFields fields =
        to_physical(state, state.velocity, compute_threads);
    GridField nonlinear_cotangent;
    for (auto& component : nonlinear_cotangent) {
        component.assign(fields.grid_cells, {});
    }
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const std::size_t cell = mode_cell(state.waves[mode], fields.grid_n);
        set_grid_value(nonlinear_cotangent, cell, output_cotangent[mode]);
    }
    transform_single_field(nonlinear_cotangent, fields.grid_n, true,
                           compute_threads);

    GridField velocity_cotangent;
    GridField vorticity_cotangent;
    for (GridField* field : {&velocity_cotangent, &vorticity_cotangent}) {
        for (auto& component : *field) {
            component.resize(fields.grid_cells);
        }
    }
    for (std::size_t cell = 0; cell < fields.grid_cells; ++cell) {
        const ComplexVector nonlinear_bar =
            grid_value(nonlinear_cotangent, cell);
        const ComplexVector velocity = grid_value(fields.velocity, cell);
        const ComplexVector vorticity = grid_value(fields.vorticity, cell);
        set_grid_value(velocity_cotangent, cell,
                       cross(nonlinear_bar, vorticity));
        set_grid_value(vorticity_cotangent, cell,
                       cross(velocity, nonlinear_bar));
    }
    transform_fields(velocity_cotangent, vorticity_cotangent,
                     fields.grid_n, false, compute_threads);

    SpectralIncrement result(state.waves.size());
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const WaveVector wave = state.waves[mode];
        const std::size_t cell = mode_cell(wave, fields.grid_n);
        const ComplexVector velocity_bar =
            grid_value(velocity_cotangent, cell);
        const ComplexVector vorticity_bar =
            grid_value(vorticity_cotangent, cell);
        const ComplexVector curl_bar =
            curl_coefficient(wave, vorticity_bar);
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                velocity_bar[component] + curl_bar[component];
        }
    }
    return result;
}

}  // namespace lemma
