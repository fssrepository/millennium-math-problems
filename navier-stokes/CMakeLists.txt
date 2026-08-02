cmake_minimum_required(VERSION 3.16)

project(navier_stokes_lemma_lab VERSION 0.1.0 LANGUAGES CXX)

option(NS_ENABLE_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)

find_package(Threads REQUIRED)
find_package(OpenMP)

add_executable(navier_stokes_lab
    src/app/main.cpp
    src/app/lemma_cli.cpp
    src/app/lemma_engine.cpp

    src/spectral/spectral_adjoint.cpp
    src/spectral/spectral_dynamics.cpp
    src/spectral/spectral_fft_operator.cpp
    src/spectral/spectral_galerkin.cpp
    src/spectral/spectral_objective.cpp
    src/spectral/spectral_state.cpp
    src/spectral/state_analysis.cpp
    src/spectral/state_transform.cpp
    src/spectral/trajectory_analyzer.cpp

    src/optimization/dynamic_adversary.cpp
    src/optimization/gradient_adversary.cpp
    src/optimization/initial_sobolev_constraint.cpp

    src/triads/local_triad_symmetrizer.cpp
    src/triads/orthogonal_triad_geometry.cpp
    src/triads/triad_commutator.cpp
    src/triads/triad_ledger.cpp
    src/triads/triad_partition.cpp
    src/triads/triad_tail_envelope.cpp
    src/triads/triad_verifier.cpp

    src/helical/helical_adversary_cli.cpp
    src/helical/helical_cutoff_scan.cpp
    src/helical/helical_gap_ledger.cpp
    src/helical/helical_sector_adjoint.cpp
    src/helical/helical_sector_adversary.cpp
    src/helical/helical_sector_objective.cpp
    src/helical/helical_trajectory_adversary.cpp
    src/helical/helical_triad_ledger.cpp

    src/proof/dyadic_shell_bounds.cpp
    src/proof/doubling_quartet_closure.cpp
    src/proof/far_tail_closure.cpp
    src/proof/moving_gap_controller.cpp
    src/proof/periodic_shell_geometry.cpp
    src/proof/periodic_tail_bound.cpp
    src/proof/proof_scaling.cpp
    src/proof/remainder_quartet_closure.cpp
    src/proof/transition_block_scaling.cpp

    src/reporting/adversary_reporter.cpp
    src/reporting/family_reporter.cpp
    src/reporting/lemma_reporter.cpp

    src/local_sld/core/local_quartic_closure_objective.cpp
    src/local_sld/core/local_signature_density.cpp
    src/local_sld/core/local_signature_objective.cpp
    src/local_sld/core/local_sld_block_objective.cpp
    src/local_sld/core/local_sld_cyclic_basis.cpp
    src/local_sld/core/local_sld_cyclic_orbit_basis.cpp
    src/local_sld/core/local_sld_response_basis.cpp
    src/local_sld/core/local_sld_two_scale_state.cpp

    src/local_sld/analysis/local_critical_derivative_ledger.cpp
    src/local_sld/analysis/local_quartic_closure_target.cpp
    src/local_sld/analysis/local_quartic_commutator.cpp
    src/local_sld/analysis/local_quartic_identity_ledger.cpp
    src/local_sld/analysis/local_quartic_projected_residual.cpp
    src/local_sld/analysis/local_quartic_reduced_ledger.cpp
    src/local_sld/analysis/local_quartic_shell_envelope.cpp
    src/local_sld/analysis/local_quartic_shell_ledger.cpp
    src/local_sld/analysis/local_signature_geometry.cpp
    src/local_sld/analysis/local_sld_response_diagonal.cpp
    src/local_sld/analysis/local_sld_doubling_shell_ledger.cpp
    src/local_sld/analysis/local_sld_doubling_scale_scan.cpp
    src/local_sld/analysis/local_sld_response_family.cpp
    src/local_sld/analysis/local_sld_response_hierarchy.cpp
    src/local_sld/analysis/local_sld_response_tensor.cpp
    src/local_sld/analysis/local_sld_projected_square.cpp
    src/local_sld/analysis/local_sld_signature_block.cpp
    src/local_sld/analysis/local_sld_trajectory_evaluator.cpp
    src/local_sld/analysis/shifted_critical_density.cpp
    src/local_sld/analysis/shifted_critical_density_budget.cpp
    src/local_sld/analysis/stretching_derivative_ledger.cpp

    src/local_sld/optimization/local_quartic_closure_adversary.cpp
    src/local_sld/optimization/local_signature_adversary.cpp
    src/local_sld/optimization/local_signature_factor_adversary.cpp
    src/local_sld/optimization/local_signature_gradient_adversary.cpp
    src/local_sld/optimization/local_signature_state_factory.cpp
    src/local_sld/optimization/local_signature_trajectory.cpp
    src/local_sld/optimization/local_sld_cyclic_ansatz.cpp
    src/local_sld/optimization/local_sld_cyclic_krylov_ansatz.cpp
    src/local_sld/optimization/local_sld_cyclic_trajectory_ansatz.cpp
    src/local_sld/optimization/local_sld_trajectory_adjoint.cpp

    src/local_sld/cli/local_quartic_closure_cli.cpp
    src/local_sld/cli/local_quartic_closure_reporter.cpp
    src/local_sld/cli/shifted_critical_density_cli.cpp
    src/local_sld/cli/shifted_critical_density_reporter.cpp
)
target_include_directories(navier_stokes_lab PRIVATE
    src/app
    src/spectral
    src/optimization
    src/triads
    src/helical
    src/proof
    src/reporting
    src/local_sld/core
    src/local_sld/analysis
    src/local_sld/optimization
    src/local_sld/cli
)
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
