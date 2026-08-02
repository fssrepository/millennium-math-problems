#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace lemma {

struct Options {
    int exponent_denominator = 32;
    int triad_cutoff = 3;
    int triad_samples = 64;
    std::uint64_t seed = 20260801;
    std::string certificate_path;
};

struct AdversaryOptions {
    std::vector<int> cutoffs{1, 2, 3};
    int restarts = 4;
    int dynamic_restarts = 1;
    int generations = 80;
    int dynamic_generations = 24;
    double mutation = 0.20;
    double viscosity = 0.10;
    double evolution_time = 0.10;
    double time_step = 0.002;
    std::uint64_t seed = 20260801;
    std::string certificate_path;
    std::string state_prefix;
    std::string state_directory;
    std::string dynamic_warm_state;
    bool dynamic_replay_each_cutoff = false;
    int sobolev_order = 0;
    double sobolev_cap = 0.0;
    double critical_density_shift = 0.0;
    std::string dynamic_objective = "critical-integral";
    std::string dynamic_optimizer = "gradient";
    std::string gradient_method = "steepest";
    int minimum_dyadic_gap = 2;
    int threads = 0;
    std::string backend = "auto";
};

struct FamilyOptions {
    std::vector<int> cutoffs{1, 2, 3, 4};
    double spectral_decay = 0.8;
    double viscosity = 0.10;
    double evolution_time = 0.10;
    double time_step = 0.002;
    std::uint64_t seed = 20260801;
    int seed_count = 1;
    std::string certificate_path;
    int threads = 0;
    std::string backend = "auto";
};

int run(const Options& options, std::ostream& out);
bool self_test(std::ostream& out);
int run_adversary(const AdversaryOptions& options, std::ostream& out);
int run_family(const FamilyOptions& options, std::ostream& out);

class LemmaCli {
public:
    static Options parse_options(int argc, char** argv, int first);
    static AdversaryOptions parse_adversary_options(int argc, char** argv,
                                                    int first);
    static FamilyOptions parse_family_options(int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static void print_adversary_help(std::ostream& out);
    static void print_family_help(std::ostream& out);
};

}  // namespace lemma
