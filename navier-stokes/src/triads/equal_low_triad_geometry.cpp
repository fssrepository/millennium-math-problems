#include "equal_low_triad_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

namespace lemma {
namespace {

bool inside_cutoff(WaveVector wave, int cutoff) {
    return std::abs(wave.x) <= cutoff &&
        std::abs(wave.y) <= cutoff &&
        std::abs(wave.z) <= cutoff;
}

}  // namespace

EqualLowTriadGeometryCertificate EqualLowTriadGeometry::certify(
    int maximum_cutoff,
    int target_length_multiplier) {
    if (maximum_cutoff < 1 || maximum_cutoff > 12) {
        throw std::invalid_argument(
            "equal-low triad cutoff must be between 1 and 12");
    }
    if (target_length_multiplier < 1 ||
        target_length_multiplier > 3) {
        throw std::invalid_argument(
            "equal-low target squared-length multiplier must be 1, 2, or 3");
    }
    EqualLowTriadGeometryCertificate certificate;
    certificate.target_length_multiplier = target_length_multiplier;
    certificate.maximum_cutoff = maximum_cutoff;
    certificate.exact_fixed_angle_relation = true;
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
        EqualLowTriadCountRow row;
        row.cutoff = cutoff;
        row.elementary_degree_bound =
            2U * static_cast<std::size_t>(2 * cutoff + 1);
        std::map<WaveVector, std::size_t> target_degree;
        for (const auto& [length, waves] : length_groups) {
            for (const WaveVector p : waves) {
                std::size_t input_degree = 0;
                for (const WaveVector q : waves) {
                    const WaveVector target = p + q;
                    if (!inside_cutoff(target, cutoff) ||
                        norm_squared(target) !=
                            target_length_multiplier * length) {
                        continue;
                    }
                    ++input_degree;
                    ++target_degree[target];
                    ++row.ordered_pairs;
                }
                row.maximum_input_degree = std::max(
                    row.maximum_input_degree, input_degree);
            }
        }
        for (const auto& [target, degree] : target_degree) {
            static_cast<void>(target);
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
