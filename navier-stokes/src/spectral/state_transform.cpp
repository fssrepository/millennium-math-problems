#include "state_transform.hpp"

#include "state_analysis.hpp"

#include <iomanip>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace lemma {

SpectralState SpectralStateTransform::to_cutoff(
    const SpectralState& state,
    int target_cutoff,
    std::uint64_t seed) {
    if (target_cutoff < 1 || target_cutoff > 12) {
        throw std::invalid_argument(
            "state transform target cutoff must be in 1..12");
    }
    const int source_cutoff = SpectralStateOps::cutoff(state);
    if (target_cutoff == source_cutoff) {
        return state;
    }
    const SpectralReal energy = SpectralStateOps::energy(state);
    SpectralState transformed;
    if (target_cutoff < source_cutoff) {
        transformed = SpectralStateFactory::project(
            state, target_cutoff);
    } else {
        std::mt19937_64 generator(seed);
        transformed = SpectralStateFactory::lift(
            state, target_cutoff, generator);
    }
    SpectralStateOps::normalize_energy(transformed, energy);
    return transformed;
}

SpectralStateTransformOptions SpectralStateTransformCli::parse(
    int argc, char** argv, int first) {
    SpectralStateTransformOptions options;
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
        } else if (name == "--output") {
            options.output_path = next(index, name);
        } else if (name == "--cutoff") {
            options.target_cutoff = std::stoi(next(index, name));
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown state-transform option: " + name);
        }
    }
    if (options.state_path.empty() || options.output_path.empty() ||
        options.target_cutoff < 1 || options.target_cutoff > 12) {
        throw std::invalid_argument(
            "state-transform requires --state, --output, and --cutoff in 1..12");
    }
    return options;
}

void SpectralStateTransformCli::print_help(std::ostream& out) {
    out << "Spectral-state cutoff transform options:\n"
        << "  --state PATH          input Fourier TSV\n"
        << "  --output PATH         output Fourier TSV\n"
        << "  --cutoff K            project or zero-pad to K\n"
        << "  --seed N              deterministic target-layout seed\n";
}

int SpectralStateTransformCli::run(
    const SpectralStateTransformOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    const int source_cutoff = SpectralStateOps::cutoff(state);
    const SpectralReal source_energy = SpectralStateOps::energy(state);
    const SpectralState transformed = SpectralStateTransform::to_cutoff(
        state, options.target_cutoff, options.seed);
    std::ostringstream metadata;
    metadata << "spectral state cutoff transform; source="
             << options.state_path
             << "; source_cutoff=" << source_cutoff
             << "; target_cutoff=" << options.target_cutoff
             << "; operation="
             << (options.target_cutoff < source_cutoff
                     ? "project"
                     : (options.target_cutoff > source_cutoff
                            ? "zero-pad"
                            : "copy"))
             << "; candidate_lemma_proved=false";
    SpectralStateWriter::write_tsv(
        options.output_path, transformed, metadata.str());
    out << std::setprecision(12)
        << "state transform " << source_cutoff << " -> "
        << options.target_cutoff
        << " energy=" << static_cast<double>(source_energy)
        << " output=" << options.output_path << '\n';
    return 0;
}

}  // namespace lemma
