cmake_minimum_required(VERSION 3.16)

project(navier_stokes_lemma_lab VERSION 0.1.0 LANGUAGES CXX)

option(NS_ENABLE_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)

find_package(Threads REQUIRED)
find_package(OpenMP)

add_executable(navier_stokes_lab
    src/dyadic_shell_bounds.cpp
    src/far_tail_closure.cpp
    src/gradient_adversary.cpp
    src/helical_triad_ledger.cpp
    src/initial_sobolev_constraint.cpp
    src/main.cpp
    src/lemma_cli.cpp
    src/lemma_engine.cpp
    src/lemma_reporter.cpp
    src/moving_gap_controller.cpp
    src/periodic_shell_geometry.cpp
    src/periodic_tail_bound.cpp
    src/adversary_reporter.cpp
    src/family_reporter.cpp
    src/proof_scaling.cpp
    src/spectral_adjoint.cpp
    src/spectral_dynamics.cpp
    src/spectral_fft_operator.cpp
    src/spectral_galerkin.cpp
    src/spectral_objective.cpp
    src/spectral_state.cpp
    src/state_analysis.cpp
    src/trajectory_analyzer.cpp
    src/transition_block_scaling.cpp
    src/triad_commutator.cpp
    src/triad_ledger.cpp
    src/triad_partition.cpp
    src/triad_tail_envelope.cpp
    src/triad_verifier.cpp
)
target_include_directories(navier_stokes_lab PRIVATE src)
target_link_libraries(navier_stokes_lab PRIVATE Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(navier_stokes_lab PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(navier_stokes_lab PRIVATE NS_HAVE_OPENMP=1)
endif()
target_compile_features(navier_stokes_lab PRIVATE cxx_std_20)
target_compile_options(navier_stokes_lab PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic;-Wconversion>
)
if(NS_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(navier_stokes_lab PRIVATE
        -fsanitize=address,undefined -fno-omit-frame-pointer
    )
    target_link_options(navier_stokes_lab PRIVATE
        -fsanitize=address,undefined
    )
endif()

include(CTest)
if(BUILD_TESTING)
    add_test(NAME navier_stokes_self_test COMMAND navier_stokes_lab self-test)
endif()
