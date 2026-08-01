#!/usr/bin/env bash
set -euo pipefail

compiler="${CXX:-g++}"
mode="${1:-release}"

flags=(-std=c++20 -Wall -Wextra -Wpedantic -Wconversion -pthread -fopenmp -DNS_HAVE_OPENMP=1)
if [[ "$mode" == "debug" ]]; then
  flags+=(-O0 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer)
elif [[ "$mode" == "release" ]]; then
  flags+=(-O3 -DNDEBUG)
else
  echo "usage: ./build.sh [release|debug]" >&2
  exit 2
fi

mkdir -p build
"$compiler" "${flags[@]}" -Isrc \
  src/gradient_adversary.cpp src/initial_sobolev_constraint.cpp \
  src/main.cpp src/lemma_cli.cpp \
  src/lemma_engine.cpp src/lemma_reporter.cpp \
  src/moving_gap_controller.cpp \
  src/adversary_reporter.cpp \
  src/family_reporter.cpp \
  src/proof_scaling.cpp src/spectral_adjoint.cpp src/spectral_dynamics.cpp \
  src/spectral_fft_operator.cpp src/spectral_galerkin.cpp \
  src/spectral_objective.cpp src/spectral_state.cpp src/state_analysis.cpp \
  src/trajectory_analyzer.cpp \
  src/triad_commutator.cpp \
  src/triad_ledger.cpp \
  src/triad_partition.cpp \
  src/triad_tail_envelope.cpp \
  src/triad_verifier.cpp \
  -o build/navier_stokes_lab
build/navier_stokes_lab self-test
