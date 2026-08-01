#pragma once

#include "local_critical_derivative_ledger.hpp"
#include "local_quartic_identity_ledger.hpp"
#include "shifted_critical_density.hpp"
#include "shifted_critical_density_budget.hpp"

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
    LocalCriticalDerivativeLedgerReport derivative_ledger;
    LocalQuarticIdentityReport quartic_identity;
    ShiftedCriticalDensityBudget derivative_budget;
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
