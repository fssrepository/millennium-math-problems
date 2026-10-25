#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
mkdir -p .verification

# This driver needs direct evaluation without the CLI optimizer's energy
# normalization. Compile the existing owners; no mathematical kernel is copied.
# Unused optimizer/FFT routines in those translation units are discarded.
"${CXX:-c++}" -std=c++20 -O1 -Wall -Wextra -Wpedantic \
  -ffunction-sections -fdata-sections -Wl,--gc-sections -pthread \
  -fopenmp -DNS_HAVE_OPENMP=1 \
  -I../src/spectral -I../src/triads -I../src/proof -I../src/proof/projective \
  -I../src/local_sld/core/closure -I../src/local_sld/core/projective \
  tests/cpp/witness.cpp \
  ../src/spectral/spectral_state.cpp \
  ../src/spectral/spectral_galerkin.cpp \
  ../src/spectral/spectral_dynamics.cpp \
  ../src/triads/triad_partition.cpp \
  ../src/triads/projective_height_shell_partition.cpp \
  ../src/triads/projective_advection_decomposition.cpp \
  ../src/proof/projective/projective_core_family.cpp \
  ../src/proof/proof_scaling.cpp \
  ../src/local_sld/core/closure/local_quartic_closure_objective.cpp \
  ../src/local_sld/core/projective/local_sld_projective_normalization_objective.cpp \
  ../src/local_sld/core/projective/local_sld_projective_normalization_cauchy_objective.cpp \
  ../src/local_sld/core/projective/local_sld_projective_normalization_schur_objective.cpp \
  -o .verification/cpp-witness
.verification/cpp-witness
