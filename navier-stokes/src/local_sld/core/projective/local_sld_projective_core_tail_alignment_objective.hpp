#pragma once

#include "projective_advection_decomposition.hpp"
#include "spectral_dynamics.hpp"

#include <cstddef>

namespace lemma {

enum class LocalSldProjectiveHeightRegion {
    core,
    tail,
};

struct LocalSldProjectiveCoreTailAlignmentObjectiveValue {
    SpectralInteger core_maximum_height = 0;
    LocalSldProjectiveHeightRegion region =
        LocalSldProjectiveHeightRegion::tail;
    std::size_t projective_shape_count = 0;
    SpectralReal signed_stretching = 0.0L;
    SpectralReal enstrophy = 0.0L;
    SpectralReal aggregate_h1_norm2 = 0.0L;
    SpectralReal stretching_h1_alignment_squared = 0.0L;
    bool finite = false;
};

class LocalSldProjectiveCoreTailAlignmentObjective {
public:
    LocalSldProjectiveCoreTailAlignmentObjective(
        const SpectralDynamics& dynamics,
        TriadSelection selection,
        SpectralInteger core_maximum_height,
        LocalSldProjectiveHeightRegion region,
        int threads = 12);

    [[nodiscard]] LocalSldProjectiveCoreTailAlignmentObjectiveValue
    evaluate(const SpectralState& state) const;

    [[nodiscard]] SpectralIncrement gradient(
        const SpectralState& state) const;

private:
    const SpectralDynamics& dynamics_;
    TriadSelection selection_;
    SpectralInteger core_maximum_height_ = 0;
    LocalSldProjectiveHeightRegion region_ =
        LocalSldProjectiveHeightRegion::tail;
    int threads_ = 12;
};

}  // namespace lemma
