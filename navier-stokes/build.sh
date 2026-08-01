#!/usr/bin/env bash
set -euo pipefail

mode="${1:-release}"

if [[ "$mode" == "debug" ]]; then
  build_type="Debug"
elif [[ "$mode" == "release" ]]; then
  build_type="Release"
else
  echo "usage: ./build.sh [release|debug]" >&2
  exit 2
fi

cmake_bin="${CMAKE:-cmake}"
if ! command -v "$cmake_bin" >/dev/null 2>&1; then
  echo "cmake not found; set CMAKE=/path/to/cmake" >&2
  exit 2
fi

build_directory="${NS_BUILD_DIR:-build/cmake-${mode}}"
build_jobs="${NS_BUILD_JOBS:-$(nproc)}"
cmake_options=(
  -DCMAKE_BUILD_TYPE="$build_type"
  -DBUILD_TESTING=ON
)
if [[ "$mode" == "debug" ]]; then
  cmake_options+=(-DNS_ENABLE_SANITIZERS=ON)
fi

"$cmake_bin" -S . -B "$build_directory" "${cmake_options[@]}"
"$cmake_bin" --build "$build_directory" --parallel "$build_jobs"
"$build_directory/navier_stokes_lab" self-test
