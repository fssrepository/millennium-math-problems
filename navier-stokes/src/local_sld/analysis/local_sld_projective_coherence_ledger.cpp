#include "local_sld_projective_coherence_ledger.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {
namespace {

using Shape = std::array<SpectralInteger, 3>;

struct ShapeInteractions {
    Shape primitive_squared_lengths{};
    std::vector<InteractionIndex> interactions;
};

Shape primitive_shape(
    const SpectralState& state,
    InteractionIndex interaction) {
    const auto [p, q, target] = interaction;
    Shape signature{
        norm_squared(state.waves[p]),
        norm_squared(state.waves[q]),
        norm_squared(state.waves[target])};
    std::sort(signature.begin(), signature.end());
    const SpectralInteger divisor = std::gcd(
        signature[0], std::gcd(signature[1], signature[2]));
    return {signature[0] / divisor,
            signature[1] / divisor,
            signature[2] / divisor};
}

TriadSelection selection_for(
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (exclude_triple_family) {
        return exclude_signature_123
            ? TriadSelection::
                  local_without_equal_low_double_triple_and_signature(
                      1, 2, 3)
            : TriadSelection::local_without_equal_low_double_triple();
    }
    return exclude_signature_123
        ? TriadSelection::
              local_without_equal_low_doubling_and_signature(1, 2, 3)
        : TriadSelection::local_without_equal_low_doubling();
}

void clear_increment(SpectralIncrement& value) {
    std::fill(value.begin(), value.end(), ComplexVector{});
}

void add_increment(
    SpectralIncrement& target,
    const SpectralIncrement& source) {
    for (std::size_t mode = 0; mode < target.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            target[mode][component] += source[mode][component];
        }
    }
}

SpectralReal increment_norm2(const SpectralIncrement& value) {
    SpectralReal result = 0.0L;
    for (const ComplexVector& mode : value) {
        for (const SpectralComplex component : mode) {
            result += std::norm(component);
        }
    }
    return result;
}

SpectralReal difference_norm2(
    const SpectralIncrement& left,
    const SpectralIncrement& right) {
    SpectralReal result = 0.0L;
    for (std::size_t mode = 0; mode < left.size(); ++mode) {
        for (std::size_t component = 0; component < 3; ++component) {
            result += std::norm(
                left[mode][component] - right[mode][component]);
        }
    }
    return result;
}

void accumulate_shape(
    const SpectralState& state,
    const ShapeInteractions& shape,
    SpectralIncrement& result) {
    clear_increment(result);
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction : shape.interactions) {
        const auto [p, q, target] = interaction;
        const SpectralComplex coefficient = imaginary_unit *
            wave_dot(state.waves[q], state.velocity[p]);
        for (std::size_t component = 0; component < 3; ++component) {
            result[target][component] +=
                coefficient * state.velocity[q][component];
        }
    }
    for (std::size_t mode = 0; mode < result.size(); ++mode) {
        result[mode] = project_divergence_free(
            state.waves[mode], result[mode]);
    }
}

bool nonzero_mode(const ComplexVector& value) {
    return std::norm(value[0]) + std::norm(value[1]) +
        std::norm(value[2]) > 1e-60L;
}

SpectralReal mode_norm2(const ComplexVector& value) {
    return std::norm(value[0]) + std::norm(value[1]) +
        std::norm(value[2]);
}

}  // namespace

LocalSldProjectiveCoherenceReport
LocalSldProjectiveCoherenceLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    int threads,
    bool exclude_signature_123,
    bool exclude_triple_family) {
    if (threads < 1 || threads > 256) {
        throw std::invalid_argument(
            "projective coherence threads must be between 1 and 256");
    }
    LocalSldProjectiveCoherenceReport report;
    report.cutoff = SpectralStateOps::cutoff(state);
    report.threads = threads;
    report.excludes_signature_123 = exclude_signature_123;
    report.excludes_triple_family = exclude_triple_family;
    const TriadSelection selection = selection_for(
        exclude_signature_123, exclude_triple_family);
    std::map<Shape, std::vector<InteractionIndex>> grouped;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (!TriadPartitioner::includes(state, interaction, selection)) {
            continue;
        }
        grouped[primitive_shape(state, interaction)].push_back(interaction);
        ++report.selected_interactions;
    }
    std::vector<ShapeInteractions> shapes;
    shapes.reserve(grouped.size());
    for (auto& [shape, interactions] : grouped) {
        shapes.push_back({shape, std::move(interactions)});
    }
    report.projective_shape_count = shapes.size();
    const SpectralIncrement total =
        dynamics.advection_direct_partition(state, selection);
    report.total_advection_norm2 = increment_norm2(total);
    const int worker_count = std::max(
        1, std::min(threads, static_cast<int>(shapes.size())));
    std::vector<SpectralIncrement> partial_totals(
        static_cast<std::size_t>(worker_count),
        SpectralIncrement(state.waves.size()));
    std::vector<std::vector<std::size_t>> partial_output_counts(
        static_cast<std::size_t>(worker_count),
        std::vector<std::size_t>(state.waves.size(), 0));
    std::vector<std::vector<SpectralReal>> partial_output_square_norm2(
        static_cast<std::size_t>(worker_count),
        std::vector<SpectralReal>(state.waves.size(), 0.0L));
    std::vector<SpectralReal> shape_norm2(shapes.size(), 0.0L);
#ifdef NS_HAVE_OPENMP
#pragma omp parallel num_threads(worker_count) if(worker_count > 1)
    {
        const int worker = omp_get_thread_num();
        SpectralIncrement component(state.waves.size());
#pragma omp for schedule(dynamic, 1)
        for (std::ptrdiff_t shape_index = 0;
             shape_index < static_cast<std::ptrdiff_t>(shapes.size());
             ++shape_index) {
            const std::size_t index = static_cast<std::size_t>(shape_index);
            accumulate_shape(state, shapes[index], component);
            shape_norm2[index] = increment_norm2(component);
            add_increment(
                partial_totals[static_cast<std::size_t>(worker)],
                component);
            auto& output_counts = partial_output_counts[
                static_cast<std::size_t>(worker)];
            auto& output_square_norm2 = partial_output_square_norm2[
                static_cast<std::size_t>(worker)];
            for (std::size_t mode = 0; mode < component.size(); ++mode) {
                output_square_norm2[mode] += mode_norm2(component[mode]);
                if (nonzero_mode(component[mode])) {
                    ++output_counts[mode];
                }
            }
        }
    }
#else
    static_cast<void>(worker_count);
    SpectralIncrement component(state.waves.size());
    for (std::size_t index = 0; index < shapes.size(); ++index) {
        accumulate_shape(state, shapes[index], component);
        shape_norm2[index] = increment_norm2(component);
        add_increment(partial_totals[0], component);
        for (std::size_t mode = 0; mode < component.size(); ++mode) {
            partial_output_square_norm2[0][mode] +=
                mode_norm2(component[mode]);
            if (nonzero_mode(component[mode])) {
                ++partial_output_counts[0][mode];
            }
        }
    }
#endif
    SpectralIncrement reconstructed(state.waves.size());
    std::vector<std::size_t> output_counts(state.waves.size(), 0);
    std::vector<SpectralReal> output_square_norm2(
        state.waves.size(), 0.0L);
    for (std::size_t worker = 0;
         worker < partial_totals.size(); ++worker) {
        add_increment(reconstructed, partial_totals[worker]);
        for (std::size_t mode = 0; mode < output_counts.size(); ++mode) {
            output_counts[mode] += partial_output_counts[worker][mode];
            output_square_norm2[mode] +=
                partial_output_square_norm2[worker][mode];
        }
    }
    report.maximum_output_shape_multiplicity = *std::max_element(
        output_counts.begin(), output_counts.end());
    report.projective_square_function_norm2 = std::accumulate(
        shape_norm2.begin(), shape_norm2.end(), 0.0L);
    SpectralReal shape_norm4_sum = 0.0L;
    SpectralReal maximum_shape_norm2 = 0.0L;
    for (const SpectralReal norm2 : shape_norm2) {
        shape_norm4_sum += norm2 * norm2;
        maximum_shape_norm2 = std::max(maximum_shape_norm2, norm2);
    }
    if (report.projective_square_function_norm2 > 0.0L) {
        report.coherent_synthesis_ratio = report.total_advection_norm2 /
            report.projective_square_function_norm2;
        report.coherent_synthesis_amplification = std::sqrt(
            report.coherent_synthesis_ratio);
        report.dominant_advection_shape_fraction =
            maximum_shape_norm2 /
            report.projective_square_function_norm2;
    }
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        if (!(output_square_norm2[mode] > 1e-60L)) {
            continue;
        }
        const SpectralReal coherent_mode_norm2 = mode_norm2(total[mode]);
        const SpectralReal ratio = coherent_mode_norm2 /
            output_square_norm2[mode];
        if (ratio > report.maximum_output_synthesis_ratio) {
            report.maximum_output_synthesis_ratio = ratio;
            report.maximum_output_synthesis_amplification =
                std::sqrt(ratio);
            report.maximum_output_wave = state.waves[mode];
            report.maximum_output_coherent_fraction =
                coherent_mode_norm2 /
                std::max(report.total_advection_norm2, 1e-60L);
            report.maximum_output_square_function_fraction =
                output_square_norm2[mode] /
                std::max(
                    report.projective_square_function_norm2, 1e-60L);
        }
    }
    if (shape_norm4_sum > 0.0L) {
        report.effective_advection_shapes =
            report.projective_square_function_norm2 *
            report.projective_square_function_norm2 / shape_norm4_sum;
    }
    report.reconstruction_relative_error = std::sqrt(
        difference_norm2(total, reconstructed) /
        std::max(report.total_advection_norm2, 1e-60L));
    report.exact_projective_reconstruction =
        report.reconstruction_relative_error < 1e-13L;
    return report;
}

LocalSldProjectiveCoherenceCliOptions
LocalSldProjectiveCoherenceCli::parse(
    int argc, char** argv, int first) {
    LocalSldProjectiveCoherenceCliOptions options;
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
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else if (name == "--exclude-123") {
            options.exclude_signature_123 = true;
        } else if (name == "--exclude-triple-family") {
            options.exclude_triple_family = true;
        } else {
            throw std::invalid_argument(
                "unknown projective-coherence option: " + name);
        }
    }
    if (options.state_path.empty() || options.certificate_path.empty() ||
        options.threads < 1 || options.threads > 256) {
        throw std::invalid_argument(
            "projective-coherence requires state, certificate, and valid threads");
    }
    return options;
}

void LocalSldProjectiveCoherenceCli::print_help(std::ostream& out) {
    out << "Local SLD projective synthesis-coherence options:\n"
        << "  --state PATH          replayable Fourier-state TSV\n"
        << "  --certificate PATH    write English JSON coherence ledger\n"
        << "  --threads N           parallel projective-ray workers\n"
        << "  --exclude-123         remove the exact (1,2,3) signature\n"
        << "  --exclude-triple-family  remove every (m,m,3m) signature\n";
}

int LocalSldProjectiveCoherenceCli::run(
    const LocalSldProjectiveCoherenceCliOptions& options,
    std::ostream& out) {
    const SpectralState state = SpectralStateReader::read_tsv(
        options.state_path);
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldProjectiveCoherenceReport report =
        LocalSldProjectiveCoherenceLedger::analyze(
            dynamics, state, options.threads,
            options.exclude_signature_123,
            options.exclude_triple_family);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write projective coherence certificate");
    }
    certificate << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-projective-coherence-v1\",\n"
        << "  \"state_path\": \"" << options.state_path << "\",\n"
        << "  \"cutoff\": " << report.cutoff << ",\n"
        << "  \"threads\": " << report.threads << ",\n"
        << "  \"selected_interactions\": "
        << report.selected_interactions << ",\n"
        << "  \"projective_shape_count\": "
        << report.projective_shape_count << ",\n"
        << "  \"maximum_output_shape_multiplicity\": "
        << report.maximum_output_shape_multiplicity << ",\n"
        << "  \"excludes_signature_123\": "
        << (report.excludes_signature_123 ? "true" : "false")
        << ",\n"
        << "  \"excludes_triple_family\": "
        << (report.excludes_triple_family ? "true" : "false")
        << ",\n"
        << "  \"total_advection_norm2\": "
        << static_cast<double>(report.total_advection_norm2) << ",\n"
        << "  \"projective_square_function_norm2\": "
        << static_cast<double>(
               report.projective_square_function_norm2) << ",\n"
        << "  \"coherent_synthesis_ratio\": "
        << static_cast<double>(report.coherent_synthesis_ratio) << ",\n"
        << "  \"coherent_synthesis_amplification\": "
        << static_cast<double>(
               report.coherent_synthesis_amplification) << ",\n"
        << "  \"maximum_output_synthesis_ratio\": "
        << static_cast<double>(
               report.maximum_output_synthesis_ratio) << ",\n"
        << "  \"maximum_output_synthesis_amplification\": "
        << static_cast<double>(
               report.maximum_output_synthesis_amplification) << ",\n"
        << "  \"maximum_output_wave\": ["
        << report.maximum_output_wave.x << ", "
        << report.maximum_output_wave.y << ", "
        << report.maximum_output_wave.z << "],\n"
        << "  \"maximum_output_coherent_fraction\": "
        << static_cast<double>(
               report.maximum_output_coherent_fraction) << ",\n"
        << "  \"maximum_output_square_function_fraction\": "
        << static_cast<double>(
               report.maximum_output_square_function_fraction) << ",\n"
        << "  \"effective_advection_shapes\": "
        << static_cast<double>(report.effective_advection_shapes)
        << ",\n"
        << "  \"dominant_advection_shape_fraction\": "
        << static_cast<double>(
               report.dominant_advection_shape_fraction) << ",\n"
        << "  \"reconstruction_relative_error\": "
        << static_cast<double>(report.reconstruction_relative_error)
        << ",\n"
        << "  \"exact_projective_reconstruction\": "
        << (report.exact_projective_reconstruction ? "true" : "false")
        << ",\n"
        << "  \"cutoff_independent_synthesis_bound_proved\": false,\n"
        << "  \"finite_ledger_is_not_a_proof\": true\n"
        << "}\n";
    out << std::setprecision(12)
        << "projective coherence cutoff=" << report.cutoff
        << " shapes=" << report.projective_shape_count
        << " effective="
        << static_cast<double>(report.effective_advection_shapes)
        << " synthesis="
        << static_cast<double>(report.coherent_synthesis_ratio)
        << " amplification="
        << static_cast<double>(report.coherent_synthesis_amplification)
        << " output_multiplicity="
        << report.maximum_output_shape_multiplicity
        << " max_output_synthesis="
        << static_cast<double>(report.maximum_output_synthesis_ratio)
        << " reconstruction="
        << static_cast<double>(report.reconstruction_relative_error)
        << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.exact_projective_reconstruction ? 0 : 2;
}

}  // namespace lemma
