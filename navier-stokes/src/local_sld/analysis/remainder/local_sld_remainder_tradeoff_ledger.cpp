#include "local_sld_remainder_tradeoff_ledger.hpp"

#include "spectral_galerkin.hpp"
#include "state_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lemma {
namespace {

SpectralReal relative_error(SpectralReal left, SpectralReal right) {
    return std::abs(left - right) /
        std::max({std::abs(left), std::abs(right), 1e-30L});
}

template <class Getter>
SpectralReal log_cutoff_slope(
    const std::vector<LocalSldRemainderTradeoffRow>& rows,
    Getter getter) {
    SpectralReal sum_x = 0.0L;
    SpectralReal sum_y = 0.0L;
    SpectralReal sum_xx = 0.0L;
    SpectralReal sum_xy = 0.0L;
    std::size_t count = 0;
    for (const LocalSldRemainderTradeoffRow& row : rows) {
        const SpectralReal value = std::abs(getter(row));
        if (row.cutoff <= 0 || !(value > 1e-30L)) {
            continue;
        }
        const SpectralReal x = std::log(
            static_cast<SpectralReal>(row.cutoff));
        const SpectralReal y = std::log(value);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        ++count;
    }
    if (count < 2) {
        return 0.0L;
    }
    const SpectralReal n = static_cast<SpectralReal>(count);
    const SpectralReal denominator = n * sum_xx - sum_x * sum_x;
    return std::abs(denominator) > 1e-30L
        ? (n * sum_xy - sum_x * sum_y) / denominator
        : 0.0L;
}

std::string json_escape(const std::string& input) {
    std::string result;
    result.reserve(input.size());
    for (const char character : input) {
        switch (character) {
            case '\\': result += "\\\\"; break;
            case '"': result += "\\\""; break;
            case '\n': result += "\\n"; break;
            case '\r': result += "\\r"; break;
            case '\t': result += "\\t"; break;
            default: result += character; break;
        }
    }
    return result;
}

void write_json(
    const LocalSldRemainderTradeoffReport& report,
    const LocalSldRemainderTradeoffCliOptions& options,
    std::ostream& output) {
    output << std::setprecision(18)
        << "{\n"
        << "  \"schema\": \"navier-stokes-local-sld-remainder-tradeoff-v1\",\n"
        << "  \"selection\": \"local signatures excluding (m,m,2m)\",\n"
        << "  \"factorization\": \"R_rem=c_rem*Phi(x), c_rem=(K_rem+G_rem)E^(1/4)/(Z^(7/4)P), x=S_full/(E^(1/4)Z^(1/4)P), Phi(x)=4x^3/(1+x^4)\",\n"
        << "  \"linear_reduction\": \"R_rem=[(K_rem+G_rem)S_full/(Z^2P^2)]*[4x^2/(1+x^4)], with 0 <= 4x^2/(1+x^4) <= 2\",\n"
        << "  \"threads\": " << options.threads << ",\n"
        << "  \"bracket_constant_cutoff_slope\": "
        << static_cast<double>(report.bracket_constant_cutoff_slope)
        << ",\n"
        << "  \"normalized_stretching_cutoff_slope\": "
        << static_cast<double>(report.normalized_stretching_cutoff_slope)
        << ",\n"
        << "  \"shape_factor_cutoff_slope\": "
        << static_cast<double>(report.shape_factor_cutoff_slope) << ",\n"
        << "  \"block_ratio_cutoff_slope\": "
        << static_cast<double>(report.block_ratio_cutoff_slope) << ",\n"
        << "  \"all_factorizations_exact\": "
        << (report.all_factorizations_exact ? "true" : "false") << ",\n"
        << "  \"cutoff_independent_tradeoff_proved\": false,\n"
        << "  \"rows\": [\n";
    for (std::size_t index = 0; index < report.rows.size(); ++index) {
        const LocalSldRemainderTradeoffRow& row = report.rows[index];
        output << "    {\n"
            << "      \"state_path\": \""
            << json_escape(options.state_paths[index]) << "\",\n"
            << "      \"cutoff\": " << row.cutoff << ",\n"
            << "      \"energy\": " << static_cast<double>(row.energy) << ",\n"
            << "      \"enstrophy\": " << static_cast<double>(row.enstrophy) << ",\n"
            << "      \"palinstrophy\": " << static_cast<double>(row.palinstrophy) << ",\n"
            << "      \"full_stretching\": " << static_cast<double>(row.full_stretching) << ",\n"
            << "      \"remainder_bracket\": " << static_cast<double>(row.remainder_bracket) << ",\n"
            << "      \"signed_lqc3_ratio\": " << static_cast<double>(row.signed_lqc3_ratio) << ",\n"
            << "      \"bracket_constant_ratio\": " << static_cast<double>(row.bracket_constant_ratio) << ",\n"
            << "      \"normalized_stretching\": " << static_cast<double>(row.normalized_stretching) << ",\n"
            << "      \"signed_linear_tradeoff\": " << static_cast<double>(row.signed_linear_tradeoff) << ",\n"
            << "      \"power_one_commutator\": " << static_cast<double>(row.power_one_commutator) << ",\n"
            << "      \"power_one_advecting\": " << static_cast<double>(row.power_one_advecting) << ",\n"
            << "      \"power_one_enstrophy_normalization\": " << static_cast<double>(row.power_one_enstrophy_normalization) << ",\n"
            << "      \"power_one_palinstrophy_normalization\": " << static_cast<double>(row.power_one_palinstrophy_normalization) << ",\n"
            << "      \"power_one_absolute_component_sum\": " << static_cast<double>(row.power_one_absolute_component_sum) << ",\n"
            << "      \"power_one_cancellation_fraction\": " << static_cast<double>(row.power_one_cancellation_fraction) << ",\n"
            << "      \"power_one_component_reconstruction_error\": " << static_cast<double>(row.power_one_component_reconstruction_error) << ",\n"
            << "      \"linear_tradeoff\": " << static_cast<double>(row.linear_tradeoff) << ",\n"
            << "      \"quadratic_tradeoff\": " << static_cast<double>(row.quadratic_tradeoff) << ",\n"
            << "      \"cubic_tradeoff\": " << static_cast<double>(row.cubic_tradeoff) << ",\n"
            << "      \"cubic_shape_numerator\": " << static_cast<double>(row.cubic_shape_numerator) << ",\n"
            << "      \"shape_denominator\": " << static_cast<double>(row.shape_denominator) << ",\n"
            << "      \"exact_shape_factor\": " << static_cast<double>(row.exact_shape_factor) << ",\n"
            << "      \"linear_shape_multiplier\": " << static_cast<double>(row.linear_shape_multiplier) << ",\n"
            << "      \"linear_reconstructed_block_ratio\": " << static_cast<double>(row.linear_reconstructed_block_ratio) << ",\n"
            << "      \"factorized_block_ratio\": " << static_cast<double>(row.factorized_block_ratio) << ",\n"
            << "      \"direct_block_ratio\": " << static_cast<double>(row.direct_block_ratio) << ",\n"
            << "      \"lqc3_to_block_depletion\": " << static_cast<double>(row.lqc3_to_block_depletion) << ",\n"
            << "      \"upper_envelope_ratio\": " << static_cast<double>(row.upper_envelope_ratio) << ",\n"
            << "      \"negative_square_ratio\": " << static_cast<double>(row.negative_square_ratio) << ",\n"
            << "      \"shape_reconstruction_error\": " << static_cast<double>(row.shape_reconstruction_error) << ",\n"
            << "      \"linear_reconstruction_error\": " << static_cast<double>(row.linear_reconstruction_error) << ",\n"
            << "      \"product_reconstruction_error\": " << static_cast<double>(row.product_reconstruction_error) << ",\n"
            << "      \"scalar_multiplier_bound_satisfied\": "
            << (row.scalar_multiplier_bound_satisfied ? "true" : "false") << ",\n"
            << "      \"exact_factorization\": "
            << (row.exact_factorization ? "true" : "false") << "\n"
            << "    }" << (index + 1 == report.rows.size() ? "\n" : ",\n");
    }
    output << "  ],\n"
        << "  \"remaining_requirement\": \"prove |(K_rem+G_rem)S_full| <= C Z^2P^2; the scalar multiplier then closes the exact remainder block with constant 2C\",\n"
        << "  \"finite_scan_is_not_a_proof\": true\n"
        << "}\n";
}

}  // namespace

LocalSldRemainderTradeoffRow
LocalSldRemainderTradeoffLedger::analyze(
    const SpectralDynamics& dynamics,
    const SpectralState& state) {
    const TriadSelection remainder =
        TriadSelection::local_without_equal_low_doubling();
    const LocalQuarticClosureObjectiveValue full =
        LocalQuarticClosureObjective(
            dynamics, TriadPartition::local).evaluate(state);
    const LocalQuarticClosureObjectiveValue selected =
        LocalQuarticClosureObjective(dynamics, remainder).evaluate(state);
    const LocalSldBlockObjectiveValue block =
        LocalSldBlockObjective(
            dynamics, remainder,
            LocalSldBlock::selected_closed).evaluate(state);
    const LocalSldRemainderDoubleSquareReport square =
        LocalSldRemainderDoubleSquare::analyze(dynamics, state);

    LocalSldRemainderTradeoffRow row;
    row.cutoff = SpectralStateOps::cutoff(state);
    row.energy = full.energy;
    row.enstrophy = full.enstrophy;
    row.palinstrophy = full.palinstrophy;
    row.full_stretching = full.signed_stretching;
    row.remainder_bracket = selected.signed_two_entry_bracket;
    row.signed_lqc3_ratio = selected.signed_lqc3_target_ratio;
    row.bracket_constant_ratio = block.block_constant_ratio;
    row.normalized_stretching = block.normalized_stretching;
    const SpectralReal absolute_x = std::abs(row.normalized_stretching);
    row.signed_linear_tradeoff = row.bracket_constant_ratio *
        row.normalized_stretching;
    const SpectralReal power_one_scale = row.full_stretching /
        (row.enstrophy * row.enstrophy *
         row.palinstrophy * row.palinstrophy);
    row.power_one_commutator = power_one_scale *
        selected.negative_commutator_pairing;
    row.power_one_advecting = power_one_scale *
        selected.advecting_slot;
    row.power_one_enstrophy_normalization = power_one_scale *
        selected.signed_stretching * selected.signed_stretching /
        (2.0L * row.enstrophy);
    row.power_one_palinstrophy_normalization = power_one_scale *
        3.0L * selected.signed_stretching *
        selected.palinstrophy_cross /
        (2.0L * row.palinstrophy);
    row.power_one_absolute_component_sum =
        std::abs(row.power_one_commutator) +
        std::abs(row.power_one_advecting) +
        std::abs(row.power_one_enstrophy_normalization) +
        std::abs(row.power_one_palinstrophy_normalization);
    if (row.power_one_absolute_component_sum > 0.0L) {
        row.power_one_cancellation_fraction = 1.0L -
            std::abs(row.signed_linear_tradeoff) /
                row.power_one_absolute_component_sum;
    }
    const SpectralReal power_one_component_sum =
        row.power_one_commutator + row.power_one_advecting +
        row.power_one_enstrophy_normalization +
        row.power_one_palinstrophy_normalization;
    row.power_one_component_reconstruction_error = relative_error(
        power_one_component_sum, row.signed_linear_tradeoff);
    row.linear_tradeoff = std::abs(row.signed_linear_tradeoff);
    row.quadratic_tradeoff = row.linear_tradeoff * absolute_x;
    row.cubic_tradeoff = row.quadratic_tradeoff * absolute_x;
    const SpectralReal x2 = row.normalized_stretching *
        row.normalized_stretching;
    const SpectralReal x3 = x2 * row.normalized_stretching;
    const SpectralReal x4 = x2 * x2;
    row.cubic_shape_numerator = 4.0L * x3;
    row.shape_denominator = 1.0L + x4;
    row.exact_shape_factor = row.cubic_shape_numerator /
        row.shape_denominator;
    row.linear_shape_multiplier = 4.0L * x2 /
        row.shape_denominator;
    row.linear_reconstructed_block_ratio =
        row.signed_linear_tradeoff * row.linear_shape_multiplier;
    row.factorized_block_ratio = row.bracket_constant_ratio *
        row.exact_shape_factor;
    row.direct_block_ratio = block.block_sld_ratio;
    if (std::abs(row.signed_lqc3_ratio) > 1e-30L) {
        row.lqc3_to_block_depletion = row.direct_block_ratio /
            row.signed_lqc3_ratio;
    }
    row.upper_envelope_ratio = square.upper_envelope_target_ratio;
    row.negative_square_ratio = square.negative_second_square_target_ratio;
    row.shape_reconstruction_error = relative_error(
        row.exact_shape_factor, block.common_shape_factor);
    row.linear_reconstruction_error = relative_error(
        row.linear_reconstructed_block_ratio,
        row.direct_block_ratio);
    row.product_reconstruction_error = relative_error(
        row.factorized_block_ratio, row.direct_block_ratio);
    row.scalar_multiplier_bound_satisfied =
        row.linear_shape_multiplier >= -1e-14L &&
        row.linear_shape_multiplier <= 2.0L + 1e-14L;
    row.exact_factorization = full.finite && selected.finite && block.finite &&
        square.exact_identity && row.shape_reconstruction_error < 1e-13L &&
        row.power_one_component_reconstruction_error < 1e-13L &&
        row.linear_reconstruction_error < 1e-13L &&
        row.product_reconstruction_error < 1e-13L &&
        row.scalar_multiplier_bound_satisfied;
    return row;
}

LocalSldRemainderTradeoffReport
LocalSldRemainderTradeoffLedger::analyze(
    const SpectralDynamics& dynamics,
    const std::vector<SpectralState>& states) {
    if (states.empty()) {
        throw std::invalid_argument(
            "remainder tradeoff ledger needs at least one state");
    }
    LocalSldRemainderTradeoffReport report;
    report.rows.reserve(states.size());
    for (const SpectralState& state : states) {
        report.rows.push_back(analyze(dynamics, state));
    }
    report.bracket_constant_cutoff_slope = log_cutoff_slope(
        report.rows, [](const auto& row) {
            return row.bracket_constant_ratio;
        });
    report.normalized_stretching_cutoff_slope = log_cutoff_slope(
        report.rows, [](const auto& row) {
            return row.normalized_stretching;
        });
    report.shape_factor_cutoff_slope = log_cutoff_slope(
        report.rows, [](const auto& row) {
            return row.exact_shape_factor;
        });
    report.block_ratio_cutoff_slope = log_cutoff_slope(
        report.rows, [](const auto& row) {
            return row.direct_block_ratio;
        });
    report.all_factorizations_exact = std::all_of(
        report.rows.begin(), report.rows.end(),
        [](const auto& row) { return row.exact_factorization; });
    return report;
}

LocalSldRemainderTradeoffCliOptions
LocalSldRemainderTradeoffCli::parse(
    int argc, char** argv, int first) {
    LocalSldRemainderTradeoffCliOptions options;
    auto next = [&](int& index, const std::string& name) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + name);
        }
        return std::string(argv[++index]);
    };
    for (int index = first; index < argc; ++index) {
        const std::string name = argv[index];
        if (name == "--state") {
            options.state_paths.push_back(next(index, name));
        } else if (name == "--certificate") {
            options.certificate_path = next(index, name);
        } else if (name == "--threads") {
            options.threads = std::stoi(next(index, name));
        } else {
            throw std::invalid_argument(
                "unknown remainder-tradeoff option: " + name);
        }
    }
    if (options.state_paths.empty() || options.certificate_path.empty() ||
        options.threads < 1) {
        throw std::invalid_argument(
            "remainder-tradeoff requires one or more states, a certificate, and positive threads");
    }
    return options;
}

void LocalSldRemainderTradeoffCli::print_help(std::ostream& out) {
    out << "Local SLD remainder bracket--shape tradeoff options:\n"
        << "  --state PATH          replayable Fourier-state TSV; repeat for a cutoff family\n"
        << "  --certificate PATH    write English JSON factor ledger\n"
        << "  --threads N           direct-kernel workers\n";
}

int LocalSldRemainderTradeoffCli::run(
    const LocalSldRemainderTradeoffCliOptions& options,
    std::ostream& out) {
    std::vector<SpectralState> states;
    states.reserve(options.state_paths.size());
    for (const std::string& path : options.state_paths) {
        states.push_back(SpectralStateReader::read_tsv(path));
    }
    SpectralGalerkin configuration;
    configuration.configure("direct", options.threads);
    const SpectralDynamics dynamics(configuration);
    const LocalSldRemainderTradeoffReport report =
        LocalSldRemainderTradeoffLedger::analyze(dynamics, states);
    const std::filesystem::path path(options.certificate_path);
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream certificate(path);
    if (!certificate) {
        throw std::runtime_error(
            "cannot write remainder tradeoff certificate");
    }
    write_json(report, options, certificate);
    out << std::setprecision(12);
    for (const LocalSldRemainderTradeoffRow& row : report.rows) {
        out << "remainder tradeoff cutoff=" << row.cutoff
            << " c=" << static_cast<double>(row.bracket_constant_ratio)
            << " x=" << static_cast<double>(row.normalized_stretching)
            << " shape=" << static_cast<double>(row.exact_shape_factor)
            << " product=" << static_cast<double>(row.direct_block_ratio)
            << " error="
            << static_cast<double>(row.product_reconstruction_error)
            << '\n';
    }
    out << "block cutoff slope="
        << static_cast<double>(report.block_ratio_cutoff_slope) << '\n'
        << "Certificate written to " << options.certificate_path << '\n';
    return report.all_factorizations_exact ? 0 : 2;
}

}  // namespace lemma
