#include "lemma_engine.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lemma {
namespace {

std::string next_value(int argc, char** argv, int& index,
                       const std::string& name) {
    if (index + 1 >= argc) {
        throw std::invalid_argument("missing value after " + name);
    }
    ++index;
    return argv[index];
}

std::vector<int> parse_cutoffs(const std::string& text,
                               const std::string& runner) {
    std::vector<int> cutoffs;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        const int cutoff = std::stoi(item);
        if (cutoff < 1 || cutoff > 12) {
            throw std::invalid_argument(runner +
                                        " cutoffs must be between 1 and 12");
        }
        cutoffs.push_back(cutoff);
    }
    std::sort(cutoffs.begin(), cutoffs.end());
    cutoffs.erase(std::unique(cutoffs.begin(), cutoffs.end()), cutoffs.end());
    if (cutoffs.empty()) {
        throw std::invalid_argument("empty " + runner + " cutoff list");
    }
    return cutoffs;
}

void validate_threads(int threads) {
    if (threads < 0 || threads > 256) {
        throw std::invalid_argument("--threads must be between 0 and 256");
    }
}

}  // namespace

Options LemmaCli::parse_options(int argc, char** argv, int first) {
    Options options;
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--exponent-denominator") {
            options.exponent_denominator =
                std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--triad-cutoff") {
            options.triad_cutoff = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--triad-samples") {
            options.triad_samples = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--seed") {
            options.seed = static_cast<std::uint64_t>(
                std::stoull(next_value(argc, argv, index, name)));
        } else if (name == "--certificate") {
            options.certificate_path = next_value(argc, argv, index, name);
        } else {
            throw std::invalid_argument("unknown lemma option: " + name);
        }
    }
    return options;
}

AdversaryOptions LemmaCli::parse_adversary_options(int argc, char** argv,
                                                   int first) {
    AdversaryOptions options;
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--cutoffs") {
            options.cutoffs =
                parse_cutoffs(next_value(argc, argv, index, name), "adversary");
        } else if (name == "--restarts") {
            options.restarts = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--generations") {
            options.generations = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--dynamic-generations") {
            options.dynamic_generations =
                std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--mutation") {
            options.mutation = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--evolve-time") {
            options.evolution_time = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--dt") {
            options.time_step = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--seed") {
            options.seed = static_cast<std::uint64_t>(
                std::stoull(next_value(argc, argv, index, name)));
        } else if (name == "--certificate") {
            options.certificate_path = next_value(argc, argv, index, name);
        } else if (name == "--state-prefix") {
            options.state_prefix = next_value(argc, argv, index, name);
        } else if (name == "--state-dir") {
            options.state_directory = next_value(argc, argv, index, name);
        } else if (name == "--dynamic-warm-state") {
            options.dynamic_warm_state = next_value(argc, argv, index, name);
        } else if (name == "--sobolev-order") {
            options.sobolev_order =
                std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--sobolev-cap") {
            options.sobolev_cap =
                std::stod(next_value(argc, argv, index, name));
        } else if (name == "--dynamic-objective") {
            options.dynamic_objective = next_value(argc, argv, index, name);
        } else if (name == "--dynamic-optimizer") {
            options.dynamic_optimizer = next_value(argc, argv, index, name);
        } else if (name == "--gradient-method") {
            options.gradient_method = next_value(argc, argv, index, name);
        } else if (name == "--minimum-dyadic-gap") {
            options.minimum_dyadic_gap =
                std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--threads") {
            options.threads = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--backend") {
            options.backend = next_value(argc, argv, index, name);
        } else {
            throw std::invalid_argument("unknown adversary option: " + name);
        }
    }
    if (options.restarts < 1 || options.restarts > 10000 ||
        options.generations < 1 || options.generations > 1000000 ||
        options.dynamic_generations < 0 ||
        options.dynamic_generations > 1000000 ||
        !(options.mutation > 0.0) || !std::isfinite(options.mutation) ||
        !(options.viscosity > 0.0) || !std::isfinite(options.viscosity) ||
        !(options.evolution_time > 0.0) ||
        !std::isfinite(options.evolution_time) ||
        !(options.time_step > 0.0) || !std::isfinite(options.time_step) ||
        options.sobolev_order < 0 || options.sobolev_order > 8 ||
        options.minimum_dyadic_gap < 1 ||
        options.minimum_dyadic_gap > 30 ||
        !(options.sobolev_cap >= 0.0) || !std::isfinite(options.sobolev_cap) ||
        ((options.sobolev_order == 0) != (options.sobolev_cap == 0.0))) {
        throw std::invalid_argument(
            "adversary numeric parameters are outside their range");
    }
    validate_threads(options.threads);
    if (options.dynamic_objective != "critical-integral" &&
        options.dynamic_objective != "critical-local-integral" &&
        options.dynamic_objective != "critical-nonlocal-integral" &&
        options.dynamic_objective != "critical-near-nonlocal-integral" &&
        options.dynamic_objective != "critical-far-nonlocal-integral" &&
        options.dynamic_objective != "critical-gap-tail-integral" &&
        options.dynamic_objective != "max-q" &&
        options.dynamic_objective != "terminal-q" &&
        options.dynamic_objective != "q-gain" &&
        options.dynamic_objective != "q-increase") {
        throw std::invalid_argument(
            "--dynamic-objective must be critical-integral, critical-local-integral, critical-nonlocal-integral, critical-near-nonlocal-integral, critical-far-nonlocal-integral, critical-gap-tail-integral, max-q, terminal-q, q-gain, or q-increase");
    }
    if (options.dynamic_optimizer != "gradient" &&
        options.dynamic_optimizer != "mutate" &&
        options.dynamic_optimizer != "hybrid") {
        throw std::invalid_argument(
            "--dynamic-optimizer must be gradient, mutate, or hybrid");
    }
    if (options.gradient_method != "steepest" &&
        options.gradient_method != "lbfgs") {
        throw std::invalid_argument(
            "--gradient-method must be steepest or lbfgs");
    }
    return options;
}

FamilyOptions LemmaCli::parse_family_options(int argc, char** argv, int first) {
    FamilyOptions options;
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--cutoffs") {
            options.cutoffs =
                parse_cutoffs(next_value(argc, argv, index, name), "family");
        } else if (name == "--decay") {
            options.spectral_decay = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--time") {
            options.evolution_time = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--dt") {
            options.time_step = std::stod(next_value(argc, argv, index, name));
        } else if (name == "--seed") {
            options.seed = static_cast<std::uint64_t>(
                std::stoull(next_value(argc, argv, index, name)));
        } else if (name == "--seed-count") {
            options.seed_count = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next_value(argc, argv, index, name);
        } else if (name == "--threads") {
            options.threads = std::stoi(next_value(argc, argv, index, name));
        } else if (name == "--backend") {
            options.backend = next_value(argc, argv, index, name);
        } else {
            throw std::invalid_argument("unknown family option: " + name);
        }
    }
    if (!(options.spectral_decay > 0.0) ||
        !std::isfinite(options.spectral_decay) ||
        !(options.viscosity > 0.0) || !std::isfinite(options.viscosity) ||
        !(options.evolution_time > 0.0) ||
        !std::isfinite(options.evolution_time) ||
        !(options.time_step > 0.0) || !std::isfinite(options.time_step)) {
        throw std::invalid_argument(
            "family numeric parameters are outside their range");
    }
    validate_threads(options.threads);
    if (options.seed_count < 1 || options.seed_count > 10000) {
        throw std::invalid_argument("--seed-count must be between 1 and 10000");
    }
    return options;
}

void LemmaCli::print_help(std::ostream& out) {
    out << "Lemma engine options:\n"
        << "  --exponent-denominator N  rational exponent scan resolution (default 32)\n"
        << "  --triad-cutoff K          Fourier cube [-K,K]^3 (default 3)\n"
        << "  --triad-samples N         random divergence-free states (default 64)\n"
        << "  --seed N                  reproducible spectral seed\n"
        << "  --certificate PATH        write machine-readable JSON certificate\n";
}

void LemmaCli::print_adversary_help(std::ostream& out) {
    out << "L4 static adversary options:\n"
        << "  --cutoffs A,B,C       Fourier cutoffs (default 1,2,3)\n"
        << "  --restarts N          independent hill climbs (default 4)\n"
        << "  --generations N       mutations per restart (default 80)\n"
        << "  --dynamic-generations N  dynamic optimizer iterations (default 24)\n"
        << "  --mutation VALUE      initial mutation radius (default 0.20)\n"
        << "  --nu VALUE            Galerkin viscosity (default 0.10)\n"
        << "  --evolve-time VALUE   dynamic L4 test horizon (default 0.10)\n"
        << "  --dt VALUE            requested RK4 time step (default 0.002)\n"
        << "  --seed N              reproducible optimizer seed\n"
        << "  --certificate PATH    write JSON result\n"
        << "  --state-prefix PATH   dump each winning Fourier state as TSV\n"
        << "  --state-dir PATH      write states under PATH/{static,dynamic}\n"
        << "  --dynamic-warm-state PATH  replay a TSV as first dynamic warm start\n"
        << "  --sobolev-order M     homogeneous initial H^M constraint\n"
        << "  --sobolev-cap VALUE   cutoff-independent squared H^M cap\n"
        << "  --dynamic-objective NAME  total/local/nonlocal/near/far/gap-tail critical integral, max-q, terminal-q, q-increase, or q-gain\n"
        << "  --minimum-dyadic-gap M  tail objective selects triads with gap >= M\n"
        << "  --dynamic-optimizer NAME  gradient, mutate, or hybrid\n"
        << "  --gradient-method NAME  steepest or projected lbfgs\n"
        << "  --threads N           worker threads; 0 uses up to 12 CPUs\n"
        << "  --backend NAME        auto, direct, or dealiased fft\n";
}

void LemmaCli::print_family_help(std::ostream& out) {
    out << "Projective smooth-family options:\n"
        << "  --cutoffs A,B,C       nested Fourier cutoffs (default 1,2,3,4)\n"
        << "  --decay VALUE         exp(-decay*|k|) coefficient decay (default 0.8)\n"
        << "  --nu VALUE            viscosity (default 0.10)\n"
        << "  --time VALUE          evolution horizon (default 0.10)\n"
        << "  --dt VALUE            RK4 time step, also halved for verification\n"
        << "  --seed N              deterministic per-wave seed\n"
        << "  --seed-count N        scan N consecutive analytic families\n"
        << "  --certificate PATH    write projective-family JSON certificate\n"
        << "  --threads N           parallel cutoff workers; 0 uses up to 12 CPUs\n"
        << "  --backend NAME        auto, direct, or dealiased fft\n";
}

}  // namespace lemma
