#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace fft3d {

template <typename Real>
void transform_1d(std::vector<std::complex<Real>>& values, bool inverse) {
    const std::size_t n = values.size();
    if (n == 0 || (n & (n - 1U)) != 0U) {
        throw std::invalid_argument("FFT length must be a nonzero power of two");
    }
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1U;
        for (; (j & bit) != 0U; bit >>= 1U) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(values[i], values[j]);
        }
    }

    const Real pi = std::acos(static_cast<Real>(-1));
    for (std::size_t length = 2; length <= n; length <<= 1U) {
        const Real angle = (inverse ? static_cast<Real>(2) : static_cast<Real>(-2)) *
                           pi / static_cast<Real>(length);
        const std::complex<Real> root{std::cos(angle), std::sin(angle)};
        for (std::size_t start = 0; start < n; start += length) {
            std::complex<Real> phase{1, 0};
            const std::size_t half = length >> 1U;
            for (std::size_t offset = 0; offset < half; ++offset) {
                const auto even = values[start + offset];
                const auto odd = values[start + offset + half] * phase;
                values[start + offset] = even + odd;
                values[start + offset + half] = even - odd;
                phase *= root;
            }
        }
    }
    if (inverse) {
        const Real inverse_n = static_cast<Real>(1) / static_cast<Real>(n);
        for (auto& value : values) {
            value *= inverse_n;
        }
    }
}

template <typename Real>
void transform_3d(std::vector<std::complex<Real>>& values, std::size_t n,
                  bool inverse) {
    if (values.size() != n * n * n) {
        throw std::invalid_argument("3D FFT storage size mismatch");
    }
    auto index = [n](std::size_t x, std::size_t y, std::size_t z) {
        return (z * n + y) * n + x;
    };
    std::vector<std::complex<Real>> line(n);

    for (std::size_t z = 0; z < n; ++z) {
        for (std::size_t y = 0; y < n; ++y) {
            for (std::size_t x = 0; x < n; ++x) {
                line[x] = values[index(x, y, z)];
            }
            transform_1d(line, inverse);
            for (std::size_t x = 0; x < n; ++x) {
                values[index(x, y, z)] = line[x];
            }
        }
    }
    for (std::size_t z = 0; z < n; ++z) {
        for (std::size_t x = 0; x < n; ++x) {
            for (std::size_t y = 0; y < n; ++y) {
                line[y] = values[index(x, y, z)];
            }
            transform_1d(line, inverse);
            for (std::size_t y = 0; y < n; ++y) {
                values[index(x, y, z)] = line[y];
            }
        }
    }
    for (std::size_t y = 0; y < n; ++y) {
        for (std::size_t x = 0; x < n; ++x) {
            for (std::size_t z = 0; z < n; ++z) {
                line[z] = values[index(x, y, z)];
            }
            transform_1d(line, inverse);
            for (std::size_t z = 0; z < n; ++z) {
                values[index(x, y, z)] = line[z];
            }
        }
    }
}

inline std::size_t next_power_of_two(std::size_t value) {
    std::size_t result = 1;
    while (result < value) {
        result <<= 1U;
    }
    return result;
}

}  // namespace fft3d
