#include "projective_core_family.hpp"

#include <numeric>
#include <stdexcept>

namespace lemma {

std::vector<ProjectivePrimitiveSignature>
ProjectiveCoreFamily::through_maximum_height(
    SpectralInteger maximum_height) {
    if (maximum_height < 1 || maximum_height > 256) {
        throw std::invalid_argument(
            "projective core maximum height must be 1..256");
    }
    std::vector<ProjectivePrimitiveSignature> result;
    for (SpectralInteger c = 1; c <= maximum_height; ++c) {
        for (SpectralInteger b = 1; b <= c; ++b) {
            for (SpectralInteger a = 1; a <= b; ++a) {
                if (std::gcd(a, std::gcd(b, c)) != 1) {
                    continue;
                }
                const SpectralInteger difference = c - a - b;
                if (difference * difference > 4 * a * b) {
                    continue;
                }
                result.push_back({a, b, c});
            }
        }
    }
    return result;
}

}  // namespace lemma
