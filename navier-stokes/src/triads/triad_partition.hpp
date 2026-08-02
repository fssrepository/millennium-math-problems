#pragma once

#include "spectral_state.hpp"

#include <array>
#include <limits>

namespace lemma {

enum class TriadPartition {
    all,
    local,
    nonlocal,
    near_nonlocal,
    far_nonlocal
};

class TriadSelection {
public:
    enum class SignatureMode {
        none,
        include,
        exclude,
        include_equal_low_doubling,
        exclude_equal_low_doubling
    };

    constexpr TriadSelection() = default;
    constexpr TriadSelection(TriadPartition partition) {
        switch (partition) {
            case TriadPartition::all:
                break;
            case TriadPartition::local:
                maximum_gap_ = 0;
                break;
            case TriadPartition::nonlocal:
                minimum_gap_ = 1;
                break;
            case TriadPartition::near_nonlocal:
                minimum_gap_ = 1;
                maximum_gap_ = 1;
                break;
            case TriadPartition::far_nonlocal:
                minimum_gap_ = 2;
                break;
        }
    }

    [[nodiscard]] static constexpr TriadSelection dyadic_tail(
        int minimum_gap) {
        return TriadSelection(minimum_gap, maximum_gap_value);
    }

    [[nodiscard]] static constexpr TriadSelection local_signature(
        SpectralInteger first, SpectralInteger second,
        SpectralInteger third) {
        return TriadSelection(
            0, 0, SignatureMode::include, {first, second, third});
    }

    [[nodiscard]] static constexpr TriadSelection local_without_signature(
        SpectralInteger first, SpectralInteger second,
        SpectralInteger third) {
        return TriadSelection(
            0, 0, SignatureMode::exclude, {first, second, third});
    }

    [[nodiscard]] static constexpr TriadSelection
    local_equal_low_doubling() {
        return TriadSelection(
            0, 0, SignatureMode::include_equal_low_doubling, {});
    }

    [[nodiscard]] static constexpr TriadSelection
    local_equal_low_doubling_shell(
        SpectralInteger minimum_low_squared,
        SpectralInteger maximum_low_squared_exclusive) {
        return TriadSelection(
            0, 0, SignatureMode::include_equal_low_doubling, {},
            minimum_low_squared, maximum_low_squared_exclusive);
    }

    [[nodiscard]] static constexpr TriadSelection
    local_without_equal_low_doubling() {
        return TriadSelection(
            0, 0, SignatureMode::exclude_equal_low_doubling, {});
    }

    [[nodiscard]] constexpr int minimum_gap() const {
        return minimum_gap_;
    }
    [[nodiscard]] constexpr int maximum_gap() const {
        return maximum_gap_;
    }
    [[nodiscard]] constexpr bool is_all() const {
        return minimum_gap_ == 0 && maximum_gap_ == maximum_gap_value &&
            signature_mode_ == SignatureMode::none;
    }
    [[nodiscard]] constexpr SignatureMode signature_mode() const {
        return signature_mode_;
    }
    [[nodiscard]] constexpr std::array<SpectralInteger, 3>
    squared_length_signature() const {
        return squared_length_signature_;
    }
    [[nodiscard]] constexpr SpectralInteger minimum_low_squared() const {
        return minimum_low_squared_;
    }
    [[nodiscard]] constexpr SpectralInteger
    maximum_low_squared_exclusive() const {
        return maximum_low_squared_exclusive_;
    }
    [[nodiscard]] constexpr bool includes_gap(int gap) const {
        return gap >= minimum_gap_ && gap <= maximum_gap_;
    }

private:
    static constexpr int maximum_gap_value =
        std::numeric_limits<int>::max();

    constexpr TriadSelection(int minimum_gap, int maximum_gap)
        : minimum_gap_(minimum_gap), maximum_gap_(maximum_gap) {}

    constexpr TriadSelection(
        int minimum_gap, int maximum_gap, SignatureMode signature_mode,
        std::array<SpectralInteger, 3> squared_length_signature,
        SpectralInteger minimum_low_squared = 0,
        SpectralInteger maximum_low_squared_exclusive =
            std::numeric_limits<SpectralInteger>::max())
        : minimum_gap_(minimum_gap), maximum_gap_(maximum_gap),
          signature_mode_(signature_mode),
          squared_length_signature_(squared_length_signature),
          minimum_low_squared_(minimum_low_squared),
          maximum_low_squared_exclusive_(
              maximum_low_squared_exclusive) {}

    int minimum_gap_ = 0;
    int maximum_gap_ = maximum_gap_value;
    SignatureMode signature_mode_ = SignatureMode::none;
    std::array<SpectralInteger, 3> squared_length_signature_{};
    SpectralInteger minimum_low_squared_ = 0;
    SpectralInteger maximum_low_squared_exclusive_ =
        std::numeric_limits<SpectralInteger>::max();
};

class TriadPartitioner {
public:
    static constexpr SpectralInteger locality_ratio_squared = 4;

    [[nodiscard]] static bool is_local(
        WaveVector first, WaveVector second, WaveVector third);
    [[nodiscard]] static int dyadic_gap(
        WaveVector first, WaveVector second, WaveVector third);
    [[nodiscard]] static bool includes(
        WaveVector first, WaveVector second, WaveVector third,
        TriadSelection selection);
    [[nodiscard]] static bool includes(
        const SpectralState& state, InteractionIndex interaction,
        TriadSelection selection);
};

}  // namespace lemma
