#include "local_quartic_closure_cli.hpp"

#include "local_quartic_closure_reporter.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>

namespace lemma {

LocalQuarticClosureAdversaryOptions LocalQuarticClosureCli::parse(
    int argc, char** argv, int first) {
    LocalQuarticClosureAdversaryOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--min-cutoff") {
            options.minimum_cutoff = std::stoi(next(index, name));
        } else if (name == "--max-cutoff") {
            options.maximum_cutoff = std::stoi(next(index, name));
        } else if (name == "--restarts") {
            options.restarts = std::stoi(next(index, name));
        } else if (name == "--workers") {
            options.workers = std::stoi(next(index, name));
        } else if (name == "--iterations") {
            options.iterations = std::stoi(next(index, name));
        } else if (name == "--line-search") {
            options.line_search_steps = std::stoi(next(index, name));
        } else if (name == "--lbfgs-history") {
            options.lbfgs_history = std::stoi(next(index, name));
        } else if (name == "--trajectory-steps") {
            options.trajectory_steps = std::stoi(next(index, name));
        } else if (name == "--nu") {
            options.viscosity = std::stold(next(index, name));
        } else if (name == "--dt") {
            options.time_step = std::stold(next(index, name));
        } else if (name == "--absorption-theta") {
            options.absorption_theta = std::stold(next(index, name));
        } else if (name == "--shape-power") {
            options.shape_power = std::stoi(next(index, name));
        } else if (name == "--step") {
            options.initial_step = std::stold(next(index, name));
        } else if (name == "--sobolev-order") {
            options.sobolev_order = std::stoi(next(index, name));
        } else if (name == "--sobolev-cap") {
            options.sobolev_cap = std::stold(next(index, name));
        } else if (name == "--seed") {
            options.seed = std::stoull(next(index, name));
        } else if (name == "--method") {
            options.method = next(index, name);
        } else if (name == "--backend") {
            options.backend = next(index, name);
        } else if (name == "--objective") {
            options.objective = next(index, name);
        } else if (name == "--selection") {
            options.selection = next(index, name);
        } else if (name == "--initial-profile") {
            options.initial_profile = next(index, name);
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--state-dir") {
            options.state_directory = next(index, name);
        } else if (name == "--warm-state") {
            options.warm_state_path = next(index, name);
        } else {
            throw std::invalid_argument(
                "unknown local-closure-adversary option: " + name);
        }
    }
    if (options.minimum_cutoff < 1 ||
        options.maximum_cutoff < options.minimum_cutoff ||
        options.maximum_cutoff > 8 || options.restarts < 1 ||
        options.workers < 1 || options.iterations < 1 ||
        options.line_search_steps < 1 ||
        options.lbfgs_history < 1 || options.lbfgs_history > 64 ||
        !(options.initial_step > 0.0L) ||
        !std::isfinite(options.initial_step) ||
        !(options.absorption_theta >= 0.0L) ||
        !(options.absorption_theta <= 1.0L) ||
        !std::isfinite(options.absorption_theta) ||
        options.shape_power < 0 || options.shape_power > 3 ||
        (options.backend != "auto" && options.backend != "direct" &&
         options.backend != "fft") ||
        (options.objective != "sld-ratio" &&
         options.objective != "closure-ratio" &&
         options.objective != "lqc3-ratio" &&
         options.objective != "signed-lqc3-ratio" &&
         options.objective != "remainder-envelope-ratio" &&
         options.objective != "remainder-absorption-ratio" &&
         options.objective != "shape-power-ratio" &&
         options.objective != "signed-closure-ratio" &&
         options.objective != "block-ratio" &&
         options.objective != "mixed-ratio" &&
         options.objective != "terminal-sld-ratio" &&
         options.objective != "maximum-sld-ratio") ||
        (options.selection != "local" &&
         options.selection != "doubling-family" &&
         options.selection != "doubling-remainder" &&
         options.selection != "remainder-without-123") ||
        (options.initial_profile != "mixed" &&
         options.initial_profile != "decaying" &&
         options.initial_profile != "flat" &&
         options.initial_profile != "outer-half-flat") ||
        ((options.objective == "block-ratio" ||
          options.objective == "mixed-ratio") &&
         (options.selection == "local" ||
          options.selection == "remainder-without-123")) ||
        ((options.objective == "terminal-sld-ratio" ||
          options.objective == "maximum-sld-ratio") &&
         (options.trajectory_steps < 1 || !(options.viscosity > 0.0L) ||
          !(options.time_step > 0.0L))) ||
        ((options.sobolev_order == 0) !=
         (options.sobolev_cap == 0.0L)) ||
        options.certificate_path.empty() ||
        options.state_directory.empty()) {
        throw std::invalid_argument(
            "local-closure-adversary requires valid search parameters, --certificate, and --state-dir");
    }
    return options;
}

void LocalQuarticClosureCli::print_help(std::ostream& out) {
    out << "Local quartic closure exact-gradient adversary options:\n"
        << "  --min-cutoff K       first Fourier cutoff\n"
        << "  --max-cutoff K       last Fourier cutoff (maximum 8)\n"
        << "  --restarts N         independent starts per cutoff\n"
        << "  --workers N          parallel restart workers (use 12)\n"
        << "  --iterations N       exact-gradient iterations per start\n"
        << "  --line-search N      backtracking trials per iteration\n"
        << "  --lbfgs-history N    limited-memory curvature pairs\n"
        << "  --trajectory-steps N RK4 steps for frozen-data trajectory objectives\n"
        << "  --nu X               viscosity for trajectory objectives\n"
        << "  --dt X               RK4 step for trajectory objectives\n"
        << "  --absorption-theta X retained first-square fraction in [0,1]\n"
        << "  --shape-power P      integer P=0..3 in |c|^2|x|^(2P)\n"
        << "  --step X             initial Riemannian step\n"
        << "  --method NAME        lbfgs or steepest\n"
        << "  --backend NAME       direct oracle, fft, or auto (default direct)\n"
        << "  --objective NAME     sld-ratio, terminal-sld-ratio, maximum-sld-ratio, lqc3-ratio, signed-lqc3-ratio, remainder-envelope-ratio, remainder-absorption-ratio, shape-power-ratio, closure-ratio, signed-closure-ratio, block-ratio, or mixed-ratio\n"
        << "  --selection NAME     local, doubling-family, doubling-remainder, or remainder-without-123\n"
        << "  --initial-profile NAME  mixed, decaying, flat, or outer-half-flat\n"
        << "  --sobolev-order M    optional homogeneous Sobolev cap\n"
        << "  --sobolev-cap X      cutoff-independent squared cap\n"
        << "  --seed N             deterministic master seed\n"
        << "  --certificate PATH   write English JSON certificate\n"
        << "  --state-dir PATH     write winning K*.tsv states\n"
        << "  --warm-state PATH    continue from a saved K or K-1 winner\n";
}

int LocalQuarticClosureCli::run(
    const LocalQuarticClosureAdversaryOptions& options,
    std::ostream& out) {
    LocalQuarticClosureAdversaryReport report =
        LocalQuarticClosureEnsemble::scan(options);
    LocalQuarticClosureReporter::write_artifacts(report, options);
    LocalQuarticClosureReporter::print_summary(report, out);
    out << "Certificate written to " << options.certificate_path
        << '\n';
    return 0;
}

}  // namespace lemma
