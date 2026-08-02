#include "spectral_state_blend.hpp"

#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lemma {
namespace {

void add_scaled(
    ComplexVector& target,
    const ComplexVector& source,
    SpectralReal factor) {
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
        target[coordinate] += factor * source[coordinate];
    }
}

void validate_state(const SpectralState& state, const char* name) {
    if (state.waves.size() != state.velocity.size() ||
        state.waves.size() != state.index.size()) {
        throw std::invalid_argument(
            std::string("state blend invalid ") + name + " layout");
    }
}

SpectralState normalized_linear_combination(
    const SpectralState& left,
    SpectralReal left_weight,
    const SpectralState& right,
    SpectralReal right_weight) {
    std::map<WaveVector, ComplexVector> values;
    for (std::size_t index = 0; index < left.waves.size(); ++index) {
        add_scaled(
            values[left.waves[index]], left.velocity[index], left_weight);
    }
    for (std::size_t index = 0; index < right.waves.size(); ++index) {
        add_scaled(
            values[right.waves[index]], right.velocity[index],
            right_weight);
    }
    SpectralState result;
    result.waves.reserve(values.size());
    result.velocity.reserve(values.size());
    for (auto& [wave, value] : values) {
        result.index.emplace(wave, result.waves.size());
        result.waves.push_back(wave);
        result.velocity.push_back(
            project_divergence_free(wave, value));
    }
    SpectralStateOps::normalize_energy(result);
    return result;
}

}  // namespace

SpectralState SpectralStateBlend::blend_on_energy_sphere(
    const SpectralState& left,
    const SpectralState& right,
    SpectralReal right_weight) {
    validate_state(left, "left");
    validate_state(right, "right");
    if (!(right_weight >= 0.0L) || !(right_weight <= 1.0L) ||
        !std::isfinite(right_weight)) {
        throw std::invalid_argument(
            "state blend right weight must be finite and in [0,1]");
    }
    const SpectralReal left_weight = std::sqrt(
        std::max(0.0L, 1.0L - right_weight * right_weight));
    return normalized_linear_combination(
        left, left_weight, right, right_weight);
}

SpectralState SpectralStateBlend::affine_normalized(
    const SpectralState& left,
    const SpectralState& right,
    SpectralReal parameter) {
    validate_state(left, "left");
    validate_state(right, "right");
    if (!std::isfinite(parameter)) {
        throw std::invalid_argument(
            "state affine parameter must be finite");
    }
    return normalized_linear_combination(
        left, 1.0L - parameter, right, parameter);
}

SpectralStateBlendOptions SpectralStateBlendCli::parse(
    int argc, char** argv, int first) {
    SpectralStateBlendOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--left-state") {
            options.left_state_path = next(index, name);
        } else if (name == "--right-state") {
            options.right_state_path = next(index, name);
        } else if (name == "--right-weight") {
            options.right_weight = std::stold(next(index, name));
        } else if (name == "--output") {
            options.output_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown state-blend option: " + name);
        }
    }
    if (options.left_state_path.empty() ||
        options.right_state_path.empty() || options.output_path.empty() ||
        !(options.right_weight >= 0.0L) ||
        !(options.right_weight <= 1.0L) ||
        !std::isfinite(options.right_weight)) {
        throw std::invalid_argument(
            "state-blend requires two states, output, and right weight in [0,1]");
    }
    return options;
}

void SpectralStateBlendCli::print_help(std::ostream& out) {
    out << "Energy-sphere spectral state blend options:\n"
        << "  --left-state PATH    primary Fourier TSV\n"
        << "  --right-state PATH   perturbing Fourier TSV\n"
        << "  --right-weight X     perturbing amplitude in [0,1]\n"
        << "  --output PATH        write normalized union-layout TSV\n";
}

int SpectralStateBlendCli::run(
    const SpectralStateBlendOptions& options,
    std::ostream& out) {
    const SpectralState left = SpectralStateReader::read_tsv(
        options.left_state_path);
    const SpectralState right = SpectralStateReader::read_tsv(
        options.right_state_path);
    const SpectralState blended =
        SpectralStateBlend::blend_on_energy_sphere(
        left, right, options.right_weight);
    const std::filesystem::path path(options.output_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ostringstream metadata;
    metadata << std::setprecision(18)
        << "energy-sphere spectral state blend; left="
        << options.left_state_path
        << "; right=" << options.right_state_path
        << "; right_weight="
        << static_cast<double>(options.right_weight)
        << "; candidate_lemma_proved=false";
    SpectralStateWriter::write_tsv(
        options.output_path, blended, metadata.str());
    out << std::setprecision(12)
        << "state blend right_weight="
        << static_cast<double>(options.right_weight)
        << " modes=" << blended.waves.size()
        << " energy="
        << static_cast<double>(SpectralStateOps::energy(blended))
        << " output=" << options.output_path << '\n';
    return 0;
}

}  // namespace lemma
