#include "projective_triad_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

bool inside_cutoff(WaveVector wave, int cutoff) {
    return std::abs(wave.x) <= cutoff &&
        std::abs(wave.y) <= cutoff &&
        std::abs(wave.z) <= cutoff;
}

using RoleLengths = std::array<SpectralInteger, 3>;

std::set<RoleLengths> ordered_roles(
    std::array<SpectralInteger, 3> signature) {
    std::sort(signature.begin(), signature.end());
    std::set<RoleLengths> roles;
    do {
        roles.insert(signature);
    } while (std::next_permutation(signature.begin(), signature.end()));
    return roles;
}

bool triangle_feasible(
    const std::array<SpectralInteger, 3>& signature) {
    const SpectralInteger a = signature[0];
    const SpectralInteger b = signature[1];
    const SpectralInteger c = signature[2];
    const SpectralInteger difference = c - a - b;
    return difference * difference <= 4 * a * b;
}

}  // namespace

ProjectiveTriadGeometryCertificate ProjectiveTriadGeometry::certify(
    int maximum_cutoff,
    std::array<SpectralInteger, 3> primitive_squared_lengths) {
    if (maximum_cutoff < 1 || maximum_cutoff > 12) {
        throw std::invalid_argument(
            "projective triad cutoff must be between 1 and 12");
    }
    std::sort(
        primitive_squared_lengths.begin(),
        primitive_squared_lengths.end());
    if (primitive_squared_lengths[0] <= 0) {
        throw std::invalid_argument(
            "projective squared lengths must be positive");
    }
    const SpectralInteger divisor = std::gcd(
        primitive_squared_lengths[0],
        std::gcd(primitive_squared_lengths[1],
                 primitive_squared_lengths[2]));
    if (divisor != 1) {
        throw std::invalid_argument(
            "projective squared lengths must be primitive");
    }
    ProjectiveTriadGeometryCertificate certificate;
    certificate.primitive_squared_lengths = primitive_squared_lengths;
    certificate.maximum_cutoff = maximum_cutoff;
    certificate.primitive_signature = true;
    certificate.triangle_feasible = triangle_feasible(
        primitive_squared_lengths);
    certificate.fixed_plane_sphere_geometry =
        certificate.triangle_feasible;
    const std::set<RoleLengths> roles = ordered_roles(
        primitive_squared_lengths);
    certificate.ordered_role_count = roles.size();

    for (int cutoff = 1; cutoff <= maximum_cutoff; ++cutoff) {
        std::map<SpectralInteger, std::vector<WaveVector>> length_groups;
        for (int z = -cutoff; z <= cutoff; ++z) {
            for (int y = -cutoff; y <= cutoff; ++y) {
                for (int x = -cutoff; x <= cutoff; ++x) {
                    if (x == 0 && y == 0 && z == 0) {
                        continue;
                    }
                    const WaveVector wave{x, y, z};
                    length_groups[norm_squared(wave)].push_back(wave);
                }
            }
        }
        ProjectiveTriadCountRow row;
        row.cutoff = cutoff;
        row.elementary_degree_bound =
            certificate.ordered_role_count * 2U *
            static_cast<std::size_t>(2 * cutoff + 1);
        std::map<WaveVector, std::size_t> input_degree;
        std::map<WaveVector, std::size_t> target_degree;
        const SpectralInteger maximum_squared =
            3 * static_cast<SpectralInteger>(cutoff) * cutoff;
        const SpectralInteger maximum_scale = maximum_squared /
            primitive_squared_lengths[2];
        for (SpectralInteger scale = 1;
             scale <= maximum_scale; ++scale) {
            for (const RoleLengths& role : roles) {
                const auto first_group = length_groups.find(
                    role[0] * scale);
                const auto second_group = length_groups.find(
                    role[1] * scale);
                if (first_group == length_groups.end() ||
                    second_group == length_groups.end()) {
                    continue;
                }
                for (const WaveVector p : first_group->second) {
                    for (const WaveVector q : second_group->second) {
                        const WaveVector target = p + q;
                        if (!inside_cutoff(target, cutoff) ||
                            norm_squared(target) != role[2] * scale) {
                            continue;
                        }
                        ++input_degree[p];
                        ++target_degree[target];
                        ++row.ordered_pairs;
                    }
                }
            }
        }
        for (const auto& [wave, degree] : input_degree) {
            static_cast<void>(wave);
            row.maximum_input_degree = std::max(
                row.maximum_input_degree, degree);
        }
        for (const auto& [wave, degree] : target_degree) {
            static_cast<void>(wave);
            row.maximum_target_degree = std::max(
                row.maximum_target_degree, degree);
        }
        certificate.maximum_input_degree_ratio = std::max(
            certificate.maximum_input_degree_ratio,
            static_cast<SpectralReal>(row.maximum_input_degree) /
                static_cast<SpectralReal>(row.elementary_degree_bound));
        certificate.maximum_target_degree_ratio = std::max(
            certificate.maximum_target_degree_ratio,
            static_cast<SpectralReal>(row.maximum_target_degree) /
                static_cast<SpectralReal>(row.elementary_degree_bound));
        certificate.all_degree_bounds_hold =
            certificate.all_degree_bounds_hold &&
            row.maximum_input_degree <= row.elementary_degree_bound &&
            row.maximum_target_degree <= row.elementary_degree_bound;
        certificate.rows.push_back(row);
    }
    return certificate;
}

}  // namespace lemma
