#pragma once

#include "shifted_critical_density.hpp"

#include <iosfwd>
#include <string>

namespace lemma {

struct ShiftedCriticalDensityOptions {
    std::string state_path;
    std::string certificate_path;
    SpectralReal viscosity = 0.1L;
    int threads = 12;
    std::string backend = "fft";
};

struct ShiftedCriticalDensityReport {
    std::string state_path;
    int cutoff = 0;
    int modes = 0;
    SpectralReal viscosity = 0.0L;
    ShiftedCriticalDensityDiagnostic diagnostic;
};

class ShiftedCriticalDensityReporter {
public:
    static void write_console(
        const ShiftedCriticalDensityReport& report, std::ostream& out);
    static void write_json(
        const ShiftedCriticalDensityReport& report, std::ostream& out);
};

class ShiftedCriticalDensityCli {
public:
    [[nodiscard]] static ShiftedCriticalDensityOptions parse(
        int argc, char** argv, int first);
    static void print_help(std::ostream& out);
    static int run(
        const ShiftedCriticalDensityOptions& options, std::ostream& out);
};

}  // namespace lemma
