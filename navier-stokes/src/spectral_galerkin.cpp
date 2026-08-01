#include "spectral_galerkin.hpp"

#include <stdexcept>

namespace lemma {

void SpectralGalerkin::configure(const std::string& backend_name,
                                 int compute_threads) {
    if (backend_name == "auto") {
        backend_ = Backend::automatic;
    } else if (backend_name == "direct") {
        backend_ = Backend::direct;
    } else if (backend_name == "fft") {
        backend_ = Backend::fft;
    } else {
        throw std::invalid_argument("backend must be auto, direct, or fft");
    }
    set_compute_threads(compute_threads);
}

void SpectralGalerkin::set_compute_threads(int compute_threads) {
    if (compute_threads < 1 || compute_threads > 256) {
        throw std::invalid_argument("compute thread count must be between 1 and 256");
    }
    compute_threads_ = compute_threads;
}

bool SpectralGalerkin::uses_fft(int cutoff) const {
    return backend_ == Backend::fft ||
           (backend_ == Backend::automatic && cutoff >= 5);
}

int SpectralGalerkin::compute_threads() const {
    return compute_threads_;
}

const char* SpectralGalerkin::backend_name() const {
    switch (backend_) {
        case Backend::automatic:
            return "auto";
        case Backend::direct:
            return "direct";
        case Backend::fft:
            return "fft";
    }
    return "unknown";
}

}  // namespace lemma
