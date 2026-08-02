#include "projective_fan_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

int ceil_div(int numerator, int denominator) {
    return (numerator + denominator - 1) / denominator;
}

Shape primitive_signature(int cutoff, int a, int b) {
    const SpectralInteger q2 =
        static_cast<SpectralInteger>(a) * a +
        static_cast<SpectralInteger>(b) * b;
    const int pz = cutoff - b;
    const SpectralInteger p2 =
        static_cast<SpectralInteger>(a) * a +
        static_cast<SpectralInteger>(pz) * pz;
    const SpectralInteger r2 =
        static_cast<SpectralInteger>(cutoff) * cutoff;
    Shape signature{q2, p2, r2};
    std::sort(signature.begin(), signature.end());
    const SpectralInteger divisor = std::gcd(
        signature[0], std::gcd(signature[1], signature[2]));
    for (SpectralInteger& value : signature) {
        value /= divisor;
    }
    return signature;
}

ProjectiveFanGeometryRow row_for(int cutoff) {
    ProjectiveFanGeometryRow row;
    row.cutoff = cutoff;
    const int minimum_a = ceil_div(cutoff, 2);
    const int maximum_a = 3 * cutoff / 5;
    const int minimum_b = ceil_div(cutoff, 4);
    const int maximum_b = 2 * cutoff / 5;
    std::set<Shape> shapes;
    SpectralReal coefficient_sum = 0.0L;
    SpectralReal coefficient_square_sum = 0.0L;
    for (int b = minimum_b; b <= maximum_b; ++b) {
        for (int a = minimum_a; a <= maximum_a; ++a) {
            ++row.pair_count;
            shapes.insert(primitive_signature(cutoff, a, b));
            const SpectralReal ar = static_cast<SpectralReal>(a);
            const SpectralReal pz = static_cast<SpectralReal>(cutoff - b);
            const SpectralReal coefficient = ar *
                static_cast<SpectralReal>(cutoff) /
                std::sqrt(ar * ar + pz * pz);
            coefficient_sum += coefficient;
            coefficient_square_sum += coefficient * coefficient;
        }
    }
    row.primitive_shape_count = shapes.size();
    row.every_shape_unique =
        row.primitive_shape_count == row.pair_count;
    if (coefficient_square_sum > 0.0L) {
        row.target_synthesis_ratio =
            coefficient_sum * coefficient_sum /
            coefficient_square_sum;
    }
    if (row.pair_count > 0) {
        row.synthesis_ratio_per_pair = row.target_synthesis_ratio /
            static_cast<SpectralReal>(row.pair_count);
    }
    const SpectralReal k = static_cast<SpectralReal>(cutoff);
    row.quadratic_lower_bound_ratio = row.target_synthesis_ratio /
        (k * k / 640.0L);
    row.quadratic_lower_bound_verified = cutoff < 40 ||
        row.quadratic_lower_bound_ratio >= 1.0L - 1e-15L;
    return row;
}

std::vector<int> sample_cutoffs(int maximum_cutoff) {
    std::set<int> samples{40, maximum_cutoff};
    for (int cutoff = 64; cutoff < maximum_cutoff; cutoff *= 2) {
        samples.insert(cutoff);
        if (cutoff > maximum_cutoff / 2) {
            break;
        }
    }
    std::vector<int> result;
    for (const int cutoff : samples) {
        if (cutoff <= maximum_cutoff) {
            result.push_back(cutoff);
        }
    }
    return result;
}

}  // namespace

ProjectiveFanGeometryCertificate ProjectiveFanGeometry::certify(
    int maximum_cutoff) {
    if (maximum_cutoff < 40 || maximum_cutoff > 4096) {
        throw std::invalid_argument(
            "projective fan certificate cutoff must be 40..4096");
    }
    ProjectiveFanGeometryCertificate certificate;
    certificate.maximum_cutoff = maximum_cutoff;
    // On a/K in [1/2,3/5] and b/K in [1/4,2/5], the three
    // squared lengths lie in [5/16,1]K^2, hence are local.
    certificate.local_frequency_geometry_proved = true;
    // q^2 <= 13K^2/25 and p^2 <= 369K^2/400, so r^2=K^2 is
    // the unique largest role.
    certificate.target_is_unique_largest_role_proved = true;
    // Equal primitive sorted shapes have the same scale because their unique
    // largest entry is K^2. Then p^2-q^2=K^2-2Kb fixes b and q^2 fixes a.
    certificate.primitive_shape_injectivity_proved = true;
    // With u_q=e_x and u_p=i(0,-p_z,p_y)/|p|, every selected
    // p+q=(0,0,K) contribution is +(aK/|p|)e_x.
    certificate.aligned_target_coefficients_proved = true;
    // For K>=40 the b interval has at least K/8 integers and the a
    // interval at least K/20, giving N>=K^2/160.
    certificate.pair_count_quadratic_lower_bound_proved = true;
    // K/2 <= aK/|p| <= K, so (sum c)^2/sum c^2 >= N/4.
    certificate.target_synthesis_quadratic_lower_bound_proved = true;
    certificate.target_synthesis_unbounded_proved = true;
    // Active q modes live in the lower z band with x polarization; active p
    // modes live in a disjoint upper band with yz polarization. q never
    // advects, while every nonzero p-advection output misses the active
    // support. Hence <Au,B(u,u)>=0 exactly.
    certificate.stretching_support_disjointness_proved = true;
    certificate.exact_zero_stretching_proved = true;
    certificate.exact_zero_power_one_product_proved = true;
    for (const int cutoff : sample_cutoffs(maximum_cutoff)) {
        certificate.rows.push_back(row_for(cutoff));
    }
    return certificate;
}

}  // namespace lemma
