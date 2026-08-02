#pragma once

#include <string>

namespace lemma {

class SpectralGalerkin {
public:
    enum class Backend { automatic, direct, fft };

    void configure(const std::string& backend_name, int compute_threads);
    void set_compute_threads(int compute_threads);

    [[nodiscard]] bool uses_fft(int cutoff) const;
    [[nodiscard]] int compute_threads() const;
    [[nodiscard]] const char* backend_name() const;

private:
    Backend backend_ = Backend::automatic;
    int compute_threads_ = 1;
};

}  // namespace lemma
