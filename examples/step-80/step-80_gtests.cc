/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function_parser.h>

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wkeyword-macro"
#endif

#define private public
#include "step-80.cc"
#undef private

#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

namespace Step80
{
  // Solve only the fluid Navier-Stokes saddle-point subproblem, ignoring the
  // immersed solid. This lives in the test rather than in the tutorial code:
  // it assumes the coupling matrices are zero (i.e. assemble_coupling() was not
  // called), so that the fluid velocity/pressure system decouples from the
  // solid displacement/multiplier blocks. It reproduces the fluid part of
  // NavierStokesImmersedProblem::solve() as the augmented Lagrangian
  // saddle-point system
  //
  //   [ A + gamma B^T Mp^{-1} B   B^T ] [u]   [f_u]
  //   [ B                          0  ] [p] = [f_p],
  //
  // solved with FGMRES preconditioned by the block-triangular AL
  // preconditioner. It accesses the problem's members directly (they are public
  // for exactly this kind of test-driven verification).
  template <int dim, int spacedim>
  void solve_navier_stokes(NavierStokesImmersedProblem<dim, spacedim> &problem,
                           const double tolerance = 0.0)
  {
    const auto &par         = problem.par;
    const auto &constraints = problem.active_fluid_constraints();

    using Vec   = LA::MPI::Vector;
    using LinOp = LinearOperator<Vec>;

    const auto A  = LinOp(problem.fluid_matrix.block(0, 0));
    const auto Bt = LinOp(problem.fluid_matrix.block(0, 1));
    const auto B  = LinOp(problem.fluid_matrix.block(1, 0));
    const auto Z6 = 0.0 * LinOp(problem.fluid_matrix.block(1, 1));
    const auto Mp = LinOp(problem.fluid_mass_matrix.block(1, 1));

    using BVec = typename LA::MPI::BlockVector;

    // Inversion of the pressure mass matrix.
    SolverControl inner_solver_control(par.inner_max_iterations,
                                       par.inner_tolerance,
                                       false,
                                       false);
    SolverCG<Vec> cg_solver(inner_solver_control);

    auto invMp =
      inverse_operator(Mp, cg_solver, problem.fluid_pressure_preconditioner);

    const auto gamma1 = par.gamma_AL_background;

    // Augmented velocity block. When the divergence penalty is already baked
    // into the assembled matrix (operator augmentation) we use A directly,
    // otherwise we add the grad-div augmentation through the pressure mass
    // matrix, exactly as in the coupled solver.
    auto A11_aug = null_operator(A);
    if (par.use_operator_augmentation)
      A11_aug = A;
    else
      A11_aug = A + gamma1 * Bt * invMp * B;

    SolverControl inner_solver_control_lagrangian(
      par.inner_lagrangian_max_iterations,
      par.inner_lagrangian_tolerance,
      false,
      par.log_inner_lagrangian_iterations);
    SolverCG<Vec> cg_solver_lagrangian(inner_solver_control_lagrangian);

    auto A11_aug_inv = inverse_operator(A11_aug,
                                        cg_solver_lagrangian,
                                        problem.fluid_velocity_preconditioner);

    std::array<std::array<LinOp, 2>, 2> system_array = {{{{A11_aug, Bt}}, // vel
                                                         {{B, Z6}}}}; // pres

    const auto system  = block_operator<2, 2, BVec>(system_array);
    auto       prec_AL = system;

    prec_AL.vmult = [&](auto &v, const auto &u) {
      v.block(0) = 0.;
      v.block(1) = 0.;

      v.block(1) = -gamma1 * invMp * u.block(1);
      v.block(0) = A11_aug_inv * (u.block(0) - Bt * v.block(1));
    };

    const std::vector<IndexSet> block_partitioning = {
      problem.fluid_owned_dofs[0], problem.fluid_owned_dofs[1]};
    BVec block_system_rhs, block_system_solution;
    block_system_rhs.reinit(block_partitioning, problem.mpi_communicator);
    block_system_rhs.block(0) = problem.fluid_system_rhs.block(0);
    block_system_rhs.block(1) = problem.fluid_system_rhs.block(1);

    block_system_solution.reinit(block_system_rhs);
    block_system_solution.block(0) =
      problem.fluid_locally_relevant_solution.block(0);
    block_system_solution.block(1) =
      problem.fluid_locally_relevant_solution.block(1);

    SolverControl                      solver_control(par.outer_max_iterations,
                                 std::max(par.outer_tolerance, tolerance),
                                 true,
                                 false);
    SolverFGMRES<LA::MPI::BlockVector> solver(solver_control);

    constraints.set_zero(problem.fluid_solution);

    solver.solve(system, block_system_solution, block_system_rhs, prec_AL);
    problem.fluid_solution.block(0) = block_system_solution.block(0);
    problem.fluid_solution.block(1) = block_system_solution.block(1);

    constraints.distribute(problem.fluid_solution);

    // Remove the constant-pressure nullspace, as in the coupled solver.
    if (problem.fluid_constant_pressure.size() == 0)
      {
        problem.fluid_constant_pressure.reinit(problem.fluid_solution);
        problem.fluid_constant_pressure.block(1) =
          invMp * problem.fluid_dual_of_constant_pressure.block(1);
      }

    const auto avg_pressure =
      problem.fluid_dual_of_constant_pressure * problem.fluid_solution;
    problem.fluid_solution.block(1).add(
      -avg_pressure, problem.fluid_constant_pressure.block(1));

    problem.fluid_locally_relevant_solution = problem.fluid_solution;
  }

  // Solve the immersed-solid elasticity saddle-point subproblem, ignoring
  // the surrounding fluid. It follows the solve_navier_stokes() above: it
  // assumes the coupling matrices are zero (i.e. assemble_coupling() was not
  // called), so that we solve for the solid displacement/multiplier system
  // as it decouples from the fluid velocity/pressure blocks. It solves the
  // augmented Lagrangian saddle-point system
  //
  //   [ K + gamma D^T M^{-1} D   D^T ] [w]   [f_w + gamma D^T M^{-1} f_l]
  //   [ D                         0  ] [l] = [f_l                      ],
  //
  // where K is the (compressible) elasticity operator, D and D^T couple the
  // displacement w to the distributed Lagrange multiplier l, and M is the
  // multiplier mass matrix used for the augmented-Lagrangian preconditioner.
  // It is solved with FGMRES preconditioned by the block-triangular AL
  // preconditioner, exactly as the solid block of the coupled solver.
  template <int dim, int spacedim>
  void solve_elasticity(NavierStokesImmersedProblem<dim, spacedim> &problem,
                        const double tolerance = 0.0)
  {
    const auto &par = problem.par;

    using Vec   = LA::MPI::Vector;
    using LinOp = LinearOperator<Vec>;

    const auto K  = LinOp(problem.solid_matrix.block(0, 0));
    const auto Dt = LinOp(problem.solid_matrix.block(0, 1));
    const auto D  = LinOp(problem.solid_matrix.block(1, 0));
    const auto Z  = 0.0 * LinOp(problem.solid_matrix.block(1, 1));
    const auto M  = LinOp(problem.solid_preconditioner.block(1, 1));

    using BVec = typename LA::MPI::BlockVector;

    // Inversion of the multiplier mass matrix.
    SolverControl inner_solver_control(par.inner_max_iterations,
                                       par.inner_tolerance,
                                       false,
                                       false);
    SolverCG<Vec> cg_solver(inner_solver_control);

    auto invW =
      inverse_operator(M, cg_solver, problem.solid_lagrange_preconditioner);

    const auto gamma2 = par.gamma_AL_immersed;

    // Augmented displacement block, mirroring the coupled solver.
    auto K_aug = K + gamma2 * Dt * invW * D;

    SolverControl inner_solver_control_lagrangian(
      par.inner_lagrangian_max_iterations,
      par.inner_lagrangian_tolerance,
      false,
      par.log_inner_lagrangian_iterations);
    SolverCG<Vec> cg_solver_lagrangian(inner_solver_control_lagrangian);

    auto K_aug_inv =
      inverse_operator(K_aug,
                       cg_solver_lagrangian,
                       problem.solid_displacement_preconditioner);

    std::array<std::array<LinOp, 2>, 2> system_array = {
      {{{K_aug, Dt}}, // displacement
       {{D, Z}}}};    // multiplier

    const auto system  = block_operator<2, 2, BVec>(system_array);
    auto       prec_AL = system;

    prec_AL.vmult = [&](auto &v, const auto &u) {
      v.block(0) = 0.;
      v.block(1) = 0.;

      v.block(1) = -gamma2 * invW * u.block(1);
      v.block(0) = K_aug_inv * (u.block(0) - Dt * v.block(1));
    };

    const std::vector<IndexSet> block_partitioning = {
      problem.solid_owned_dofs[0], problem.solid_owned_dofs[1]};
    BVec block_system_rhs, block_system_solution;
    block_system_rhs.reinit(block_partitioning, problem.mpi_communicator);
    block_system_rhs.block(0) = problem.solid_system_rhs.block(0);
    block_system_rhs.block(1) = problem.solid_system_rhs.block(1);

    // Augment the displacement right-hand side consistently with K_aug, so that
    // the augmented system is equivalent to the original saddle point.
    block_system_rhs.block(0) += gamma2 * Dt * invW * block_system_rhs.block(1);

    block_system_solution.reinit(block_system_rhs);
    block_system_solution.block(0) =
      problem.solid_locally_relevant_solution.block(0);
    block_system_solution.block(1) =
      problem.solid_locally_relevant_solution.block(1);

    SolverControl                      solver_control(par.outer_max_iterations,
                                 std::max(par.outer_tolerance, tolerance),
                                 true,
                                 false);
    SolverFGMRES<LA::MPI::BlockVector> solver(solver_control);

    solver.solve(system, block_system_solution, block_system_rhs, prec_AL);
    problem.solid_solution.block(0) = block_system_solution.block(0);
    problem.solid_solution.block(1) = block_system_solution.block(1);

    problem.solid_constraints.distribute(problem.solid_solution);
    problem.solid_locally_relevant_solution = problem.solid_solution;
  }
} // namespace Step80

namespace
{
  struct CouplingProjectionConfig
  {
    std::string id;
    std::string velocity_fe;
    std::string pressure_fe;
    std::string displacement_fe;
    std::string lagrange_fe;
    double      tolerance = 1e-11;
  };

  void run_coupling_projection_test(const CouplingProjectionConfig &config)
  {
    using namespace dealii;
    using namespace Step80;

    ParameterAcceptor::clear();

    constexpr int dim      = 2;
    constexpr int spacedim = 2;

    const std::string velocity_expr = "1 + x + 2*y; -0.5 + 3*x - y; 0";
    const std::string output_dir =
      "test-output-coupling-projection-" + config.id;

    std::ostringstream prm;
    prm << "subsection Navier-Stokes Immersed Problem\n"
        << "  set Output directory     = " << output_dir << "\n"
        << "  set Output frequency     = 1\n"
        << "  subsection Finite element spaces\n"
        << "    set Velocity            = " << config.velocity_fe << "\n"
        << "    set Pressure            = " << config.pressure_fe << "\n"
        << "    set Displacement        = " << config.displacement_fe << "\n"
        << "    set Lagrange multiplier = " << config.lagrange_fe << "\n"
        << "  end\n"
        << "  subsection Grid generation\n"
        << "    set Fluid grid generator           = hyper_cube\n"
        << "    set Fluid grid generator arguments = 0: 1: false\n"
        << "    set Initial fluid refinement       = 3\n"
        << "    set Solid grid generator           = hyper_cube\n"
        << "    set Solid grid generator arguments = 0.1: 0.9: false\n"
        << "    set Initial solid refinement       = 3\n"
        << "  end\n"
        << "  subsection Navier-Stokes boundary conditions\n"
        << "    set Function expression = " << velocity_expr << "\n"
        << "  end\n"
        << "  subsection Navier-Stokes initial conditions\n"
        << "    set Function expression = " << velocity_expr << "\n"
        << "  end\n"
        << "  subsection Solver parameters\n"
        << "    set Inner solver tolerance          = 1e-14\n"
        << "    set Inner solver maximum iterations = 2000\n"
        << "    set Coupling quadrature iterations  = 5\n"
        << "    set Coupling quadrature type        = gauss_lobatto\n"
        << "  end\n"
        << "end\n";

    const std::string prm_filename =
      (std::filesystem::temp_directory_path() /
       ("step-80-coupling-projection-" + config.id + ".prm"))
        .string();
    {
      std::ofstream out(prm_filename);
      out << prm.str();
    }

    const auto [file_dim, file_spacedim] =
      get_dimension_and_spacedimension(prm_filename);
    ASSERT_EQ(file_dim, static_cast<unsigned int>(dim));
    ASSERT_EQ(file_spacedim, static_cast<unsigned int>(spacedim));

    NavierStokesImmersedProblemParameters<dim, spacedim> par;
    ParameterAcceptor::initialize(prm_filename);

    NavierStokesImmersedProblem<dim, spacedim> problem(par);
    problem.make_grid();
    problem.initial_setup();
    problem.setup_dofs();
    problem.interpolate_initial_conditions();
    problem.setup_coupling();
    problem.assemble_coupling();
    problem.assemble_elasticity_system(1.0);

    LA::MPI::Vector projected_rhs(problem.solid_owned_dofs[1],
                                  problem.mpi_communicator);
    LA::MPI::Vector projected_multiplier(problem.solid_owned_dofs[1],
                                         problem.mpi_communicator);

    problem.coupling_matrix.block(1, 0).vmult(projected_rhs,
                                              problem.fluid_solution.block(0));

    SolverControl             solver_control(par.inner_max_iterations,
                                 par.inner_tolerance,
                                 false,
                                 false);
    SolverCG<LA::MPI::Vector> cg(solver_control);
    cg.solve(problem.solid_preconditioner.block(1, 1),
             projected_multiplier,
             projected_rhs,
             problem.solid_lagrange_preconditioner);

    problem.solid_solution                  = 0;
    problem.solid_solution.block(1)         = projected_multiplier;
    problem.solid_locally_relevant_solution = problem.solid_solution;

    problem.output_results(/*cycle=*/0, /*time=*/0.0);

    FunctionParser<spacedim> exact_multiplier(2 * spacedim);
    exact_multiplier.initialize(
      FunctionParser<spacedim>::default_variable_names(),
      "0; 0; 1 + x + 2*y; -0.5 + 3*x - y",
      {});

    const ComponentSelectFunction<spacedim> multiplier_mask(
      std::pair<unsigned int, unsigned int>(spacedim, 2 * spacedim),
      2 * spacedim);

    Vector<double> cellwise_error(
      problem.solid_dh.get_triangulation().n_active_cells());
    VectorTools::integrate_difference(problem.solid_dh,
                                      problem.solid_locally_relevant_solution,
                                      exact_multiplier,
                                      cellwise_error,
                                      QGauss<spacedim>(
                                        problem.solid_fe->degree + 2),
                                      VectorTools::L2_norm,
                                      &multiplier_mask);

    const double l2_error =
      VectorTools::compute_global_error(problem.solid_dh.get_triangulation(),
                                        cellwise_error,
                                        VectorTools::L2_norm);

    EXPECT_LT(l2_error, config.tolerance)
      << "Configuration " << config.id << " failed with FE tuple ("
      << config.velocity_fe << ", " << config.displacement_fe << ", "
      << config.lagrange_fe << ") and L2 error " << l2_error;
  }
} // namespace

TEST(Step80, GetDimensionAndSpacedimensionReadsParameterHandler)
{
  dealii::ParameterHandler prm;

  prm.declare_entry("dimension", "2", dealii::Patterns::Integer(2, 3));
  prm.declare_entry("space dimension", "2", dealii::Patterns::Integer(2, 3));
  prm.set("dimension", "3");
  prm.set("space dimension", "3");

  const auto [dim, spacedim] = Step80::get_dimension_and_spacedimension(prm);

  EXPECT_EQ(dim, 3U);
  EXPECT_EQ(spacedim, 3U);
}

TEST(Step80, MPIIsInitialized)
{
  int initialized = 0;
  ASSERT_EQ(MPI_Initialized(&initialized), MPI_SUCCESS);
  EXPECT_NE(initialized, 0);

  EXPECT_EQ(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 1U);
}

// This test exercises the actual step-80 coupling operators. We interpolate a
// manufactured velocity field into the fluid velocity space, apply the
// velocity-to-Lagrange-multiplier coupling matrix to obtain the corresponding
// right-hand side in the multiplier space, invert the multiplier mass matrix,
// and then measure the L2 error of the recovered multiplier field against the
// original analytical function on the solid mesh. We also write the fluid and
// solid fields through NavierStokesImmersedProblem::output_results() so they
// can be inspected in ParaView with the standard step-80 output layout.
class CouplingProjectionParameterizedTest
  : public ::testing::TestWithParam<CouplingProjectionConfig>
{};

TEST_P(CouplingProjectionParameterizedTest, RecoversManufacturedVelocity)
{
  run_coupling_projection_test(GetParam());
}

INSTANTIATE_TEST_SUITE_P(
  FECombinations,
  CouplingProjectionParameterizedTest,
  ::testing::Values(CouplingProjectionConfig{"q1_q1_q1",
                                             "FE_Q<2>(1)",
                                             "FE_DGP<2>(0)",
                                             "FE_Q<2>(1)",
                                             "FE_Q<2>(1)"},
                    CouplingProjectionConfig{"q2_q2_q2",
                                             "FE_Q<2>(2)",
                                             "FE_DGP<2>(1)",
                                             "FE_Q<2>(2)",
                                             "FE_Q<2>(2)"},
                    CouplingProjectionConfig{"q2_q1_q1",
                                             "FE_Q<2>(2)",
                                             "FE_DGP<2>(1)",
                                             "FE_Q<2>(1)",
                                             "FE_Q<2>(1)"},
                    CouplingProjectionConfig{"q2_q2_q1",
                                             "FE_Q<2>(2)",
                                             "FE_DGP<2>(1)",
                                             "FE_Q<2>(2)",
                                             "FE_Q<2>(1)"}),
  [](const testing::TestParamInfo<CouplingProjectionConfig> &info) {
    return info.param.id;
  });

// Method of manufactured solutions (MMS) verification of the Navier-Stokes
// part of the solver without any solid. We solve a
// steady problem on the fluid domain only, using solve_navier_stokes() (which
// ignores the solid/coupling blocks), and check that the L2 velocity error
// converges at the expected rate for the chosen velocity element (FE_Q(2)).
//
// The manufactured, steady, incompressible solution on [-1, 1]^2 with
// density = viscosity = 1 is
//
//   u = ( sin(pi x) cos(pi y), -cos(pi x) sin(pi y) )   (divergence free),
//   p =   sin(pi x) sin(pi y)                            (zero mean).
//
//
//   f_x = (pi/2) sin(2 pi x) + 2 pi^2 sin(pi x) cos(pi y) + pi cos(pi x) sin(pi
//   y), f_y = (pi/2) sin(2 pi y) - 2 pi^2 cos(pi x) sin(pi y) + pi sin(pi x)
//   cos(pi y).
TEST(Step80, NavierStokesManufacturedSolutionConvergence)
{
  using namespace dealii;
  using namespace Step80;

  constexpr int dim      = 2;
  constexpr int spacedim = 2;

  const std::string exact_expr =
    "sin(pi*x)*cos(pi*y); -cos(pi*x)*sin(pi*y); sin(pi*x)*sin(pi*y)";
  const std::string rhs_expr =
    "0.5*pi*sin(2*pi*x) + 2*pi*pi*sin(pi*x)*cos(pi*y) + "
    "pi*cos(pi*x)*sin(pi*y);"
    " 0.5*pi*sin(2*pi*y) - 2*pi*pi*cos(pi*x)*sin(pi*y) + "
    "pi*sin(pi*x)*cos(pi*y);"
    " 0";

  // Assemble a parameter file. Everything not listed here keeps its default;
  // in particular the velocity element stays FE_Q<2>(2). Solver tolerances are
  // are very low to measure convergence.
  // In all honesty, that prm generation part was vibe coded.
  std::ostringstream prm;
  prm << "subsection Navier-Stokes Immersed Problem\n"
      << "  set Dirichlet boundary ids = 0\n"
      << "  subsection Grid generation\n"
      << "    set Fluid grid generator           = hyper_cube\n"
      << "    set Fluid grid generator arguments = -1: 1: false\n"
      << "    set Initial fluid refinement       = 2\n"
      << "  end\n"
      << "  subsection Physical properties\n"
      << "    set Density   = 1\n"
      << "    set Viscosity = 1\n"
      << "  end\n"
      << "  subsection Navier-Stokes boundary conditions\n"
      << "    set Function expression = " << exact_expr << "\n"
      << "  end\n"
      << "  subsection Navier-Stokes right hand side\n"
      << "    set Function expression = " << rhs_expr << "\n"
      << "  end\n"
      << "  subsection Navier-Stokes analytical solution\n"
      << "    set Function expression = " << exact_expr << "\n"
      << "  end\n"
      << "  subsection Navier-Stokes initial conditions\n"
      << "    set Function expression = 0; 0; 0\n"
      << "  end\n"
      << "  subsection Solver parameters\n"
      << "    set Inner solver tolerance                     = 1e-12\n"
      << "    set Inner solver maximum iterations            = 1000\n"
      << "    set Outer solver tolerance                     = 1e-10\n"
      << "    set Outer solver maximum iterations            = 2000\n"
      << "    set Inner Lagrangian solver tolerance          = 1e-11\n"
      << "    set Inner Lagrangian solver maximum iterations = 3000\n"
      << "    set Gamma AL background                        = 100\n"
      << "  end\n"
      << "end\n";

  // Write the parameters to a temporary file and initialize
  const std::string prm_filename =
    (std::filesystem::temp_directory_path() / "step-80-mms.prm").string();
  {
    std::ofstream out(prm_filename);
    out << prm.str();
  }

  const auto [file_dim, file_spacedim] =
    get_dimension_and_spacedimension(prm_filename);
  ASSERT_EQ(file_dim, static_cast<unsigned int>(dim));
  ASSERT_EQ(file_spacedim, static_cast<unsigned int>(spacedim));

  NavierStokesImmersedProblemParameters<dim, spacedim> par;
  ParameterAcceptor::initialize(prm_filename);

  const std::vector<unsigned int> levels = {2, 3, 4};
  std::vector<double>             h_values;
  std::vector<double>             l2_errors;

  ConvergenceTable convergence_table;

  for (const unsigned int level : levels)
    {
      par.initial_fluid_refinement = level;

      NavierStokesImmersedProblem<dim, spacedim> problem(par);
      problem.make_grid();
      problem.initial_setup();
      problem.setup_dofs();
      problem.interpolate_initial_conditions();

      // Steady state via Picard iteration. With alpha = 0 the time term drops
      // out, so the system matrix is independent of the previous iterate and is
      // assembled once; only the (lagged) convective term in the right-hand
      // side changes from one iteration to the next. Iterate until the velocity
      // update is well below the discretization error.
      problem.assemble_navier_stokes_system(0.0);

      double change = 1.0;
      for (unsigned int it = 0; it < 50 && change > 1e-9; ++it)
        {
          LA::MPI::Vector previous_velocity = problem.fluid_solution.block(0);

          problem.assemble_navier_stokes_rhs(0.0);
          solve_navier_stokes(problem);
          problem.fluid_locally_relevant_solution_old =
            problem.fluid_locally_relevant_solution;

          previous_velocity -= problem.fluid_solution.block(0);
          change = previous_velocity.l2_norm();
        }

      // L2 error of the velocity components only (mask out the pressure).
      const ComponentSelectFunction<spacedim> velocity_mask(
        std::pair<unsigned int, unsigned int>(0, spacedim), spacedim + 1);

      Vector<double> cellwise_error(
        problem.fluid_dh.get_triangulation().n_active_cells());
      VectorTools::integrate_difference(problem.fluid_dh,
                                        problem.fluid_locally_relevant_solution,
                                        par.navier_stokes_analytical_solution,
                                        cellwise_error,
                                        QGauss<spacedim>(
                                          problem.fluid_fe->degree + 2),
                                        VectorTools::L2_norm,
                                        &velocity_mask);

      const double l2_error =
        VectorTools::compute_global_error(problem.fluid_dh.get_triangulation(),
                                          cellwise_error,
                                          VectorTools::L2_norm);

      const double h = 1.0 / std::pow(2.0, level);
      h_values.push_back(h);
      l2_errors.push_back(l2_error);

      convergence_table.add_value("level", level);
      convergence_table.add_value(
        "cells", problem.fluid_dh.get_triangulation().n_active_cells());
      convergence_table.add_value("dofs", problem.fluid_dh.n_dofs());
      convergence_table.add_value("h", h);
      convergence_table.add_value("L2_velocity", l2_error);
    }

  // Format and print the convergence table. The rate column is computed
  // assuming h halves between consecutive rows (one global refinement each).
  convergence_table.set_precision("h", 4);
  convergence_table.set_scientific("h", true);
  convergence_table.set_precision("L2_velocity", 6);
  convergence_table.set_scientific("L2_velocity", true);
  convergence_table.evaluate_convergence_rates(
    "L2_velocity", ConvergenceTable::reduction_rate_log2);

  std::cout << "\nNavier-Stokes manufactured-solution convergence:\n";
  convergence_table.write_text(std::cout);
  std::cout << std::endl;

  // For FE_Q(2) velocity we expect an L2 convergence rate of ~3. Require the
  // observed rate on each successive mesh pair to exceed velocity_degree + 0.5.
  for (unsigned int i = 1; i < levels.size(); ++i)
    {
      const double rate = std::log(l2_errors[i - 1] / l2_errors[i]) /
                          std::log(h_values[i - 1] / h_values[i]);
      EXPECT_GT(rate, 2.5)
        << "L2 velocity convergence rate between refinement levels "
        << levels[i - 1] << " and " << levels[i] << " is only " << rate
        << " (errors " << l2_errors[i - 1] << " -> " << l2_errors[i] << ").";
    }
}

// Method of manufactured solutions (MMS) verification of the solid elasticity
// part of the solver without any fluid. We solve the immersed-solid
// saddle-point problem on the solid domain only, using solve_elasticity() above
// (ignoring the fluid/coupling blocks), and check that the L2 displacement
// error converges at the expected rate for Q_1 elements.
//
// The solid subproblem couples the displacement w to a distributed Lagrange
// multiplier l through the two equations (with alpha = 1/dt, w_old = 0 here)
//
//   a(w, v) - (l, v) = (f_w, v)      for all v      [compressible elasticity]
//   -alpha (mu, w)   = (f_l, mu)     for all mu     [multiplier constraint]
//
// where a(w, v) = integral(2 mu_lame eps(w):eps(v) + lambda_lame div w div v).

// MMS: we manufacture a smooth displacement that vanishes on the solid boundary
// (so that homogeneous Dirichlet conditions make the natural boundary term drop
// out of the weak form), together with a trivial multiplier l = 0. On [-1, 1]^2
// with mu_lame = lambda_lame = 1,
//
//   w = ( sin(pi x) sin(pi y), sin(2 pi x) sin(pi y) ),   l = 0.
//
// The elasticity forcing f_w = -div sigma(w) is then:
//
//   f_w,x = 4 pi^2 sin(pi x) sin(pi y) - 4 pi^2 cos(2 pi x) cos(pi y),
//   f_w,y = 7 pi^2 sin(2 pi x) sin(pi y) - 2 pi^2 cos(pi x) cos(pi y),
//
// and the multiplier forcing f_l = -alpha w (with alpha = 1) makes the
// constraint equation reproduce (mu, w) = (mu, w_exact).
TEST(Step80, ElasticityManufacturedSolutionConvergence)
{
  using namespace dealii;
  using namespace Step80;

  // need to clear, otherwise the parameters from the previous test are still
  // registered
  ParameterAcceptor::clear();

  constexpr int dim      = 2;
  constexpr int spacedim = 2;

  // Displacement forcing f_w (components 0...spacedim-1) and multiplier forcing
  // f_l (components spacedim...2*spacedim-1), assembled into the solid right
  // hand side function expected by par.solid_rhs (2*spacedim components).
  const std::string rhs_expr =
    "4*pi*pi*sin(pi*x)*sin(pi*y) - 4*pi*pi*cos(2*pi*x)*cos(pi*y);"
    " 7*pi*pi*sin(2*pi*x)*sin(pi*y) - 2*pi*pi*cos(pi*x)*cos(pi*y);"
    " -sin(pi*x)*sin(pi*y);"
    " -sin(2*pi*x)*sin(pi*y)";

  // Assemble a parameter file. Everything not listed here keeps its default; in
  // particular the displacement element stays FE_Q<2>(1) and the multiplier
  // stays FE_DGQ<2>(0). We have a coarse fluid mesh here, since the fluid
  // system is set up by setup_dofs() but never solved in the test.
  std::ostringstream prm;
  prm << "subsection Navier-Stokes Immersed Problem\n"
      << "  subsection Grid generation\n"
      << "    set Fluid grid generator           = hyper_cube\n"
      << "    set Fluid grid generator arguments = -1: 1: false\n"
      << "    set Initial fluid refinement       = 1\n"
      << "    set Solid grid generator           = hyper_cube\n"
      << "    set Solid grid generator arguments = -1: 1: false\n"
      << "    set Initial solid refinement       = 3\n"
      << "  end\n"
      << "  subsection Physical properties\n"
      << "    set Density     = 1\n"
      << "    set Viscosity   = 1\n"
      << "    set Lame mu     = 1\n"
      << "    set Lame lambda = 1\n"
      << "  end\n"
      << "  subsection Solid right hand side\n"
      << "    set Function expression = " << rhs_expr << "\n"
      << "  end\n"
      << "  subsection Solver parameters\n"
      << "    set Inner solver tolerance                     = 1e-12\n"
      << "    set Inner solver maximum iterations            = 1000\n"
      << "    set Outer solver tolerance                     = 1e-10\n"
      << "    set Outer solver maximum iterations            = 2000\n"
      << "    set Inner Lagrangian solver tolerance          = 1e-11\n"
      << "    set Inner Lagrangian solver maximum iterations = 3000\n"
      << "    set Gamma AL immersed                          = 10\n"
      << "  end\n"
      << "end\n";

  // Write the parameters to a temporary file and initialize.
  const std::string prm_filename =
    (std::filesystem::temp_directory_path() / "step-80-solid-mms.prm").string();
  {
    std::ofstream out(prm_filename);
    out << prm.str();
  }

  const auto [file_dim, file_spacedim] =
    get_dimension_and_spacedimension(prm_filename);
  ASSERT_EQ(file_dim, static_cast<unsigned int>(dim));
  ASSERT_EQ(file_spacedim, static_cast<unsigned int>(spacedim));

  NavierStokesImmersedProblemParameters<dim, spacedim> par;
  ParameterAcceptor::initialize(prm_filename);

  // Exact displacement (plus a zero multiplier) used for the error evaluation.
  FunctionParser<spacedim> exact_solution(2 * spacedim);
  exact_solution.initialize(FunctionParser<spacedim>::default_variable_names(),
                            "sin(pi*x)*sin(pi*y); sin(2*pi*x)*sin(pi*y); 0; 0",
                            {{"pi", numbers::PI}});

  const std::vector<unsigned int> levels = {3, 4, 5};
  std::vector<double>             h_values;
  std::vector<double>             l2_errors;

  ConvergenceTable convergence_table;

  for (const unsigned int level : levels)
    {
      par.initial_solid_refinement = level;

      NavierStokesImmersedProblem<dim, spacedim> problem(par);
      problem.make_grid();
      problem.initial_setup();
      problem.setup_dofs();

      // Impose homogeneous Dirichlet conditions on the solid displacement (the
      // manufactured displacement vanishes on the boundary).
      // These conditions are consistent and make the natural boundary term
      // drop out of the weak form. setup_dofs() leaves solid_constraints
      // holding hanging-node constraints only, so we rebuild it here.
      {
        const IndexSet locally_relevant_solid_dofs =
          DoFTools::extract_locally_relevant_dofs(problem.solid_dh);
        problem.solid_constraints.clear();
        problem.solid_constraints.reinit(locally_relevant_solid_dofs,
                                         locally_relevant_solid_dofs);
        DoFTools::make_hanging_node_constraints(problem.solid_dh,
                                                problem.solid_constraints);
        VectorTools::interpolate_boundary_values(
          problem.solid_dh,
          0,
          Functions::ZeroFunction<spacedim>(2 * spacedim),
          problem.solid_constraints,
          problem.solid_fe->component_mask(problem.displacement));
        problem.solid_constraints.close();

        // Rebuild the solid matrix and preconditioner sparsity patterns so that
        // they reserve the diagonal entries needed by the (now larger) set of
        // constrained degrees of freedom, exactly as setup_dofs() does.
        const auto locally_owned_solid_dofs_per_processor =
          Utilities::MPI::all_gather(problem.mpi_communicator,
                                     problem.solid_dh.locally_owned_dofs());
        {
          BlockDynamicSparsityPattern dsp(problem.solid_dofs_per_block,
                                          problem.solid_dofs_per_block);
          DoFTools::make_sparsity_pattern(problem.solid_dh,
                                          dsp,
                                          problem.solid_constraints,
                                          false);
          SparsityTools::distribute_sparsity_pattern(
            dsp,
            locally_owned_solid_dofs_per_processor,
            problem.mpi_communicator,
            locally_relevant_solid_dofs);
          problem.solid_matrix.reinit(problem.solid_owned_dofs,
                                      dsp,
                                      problem.mpi_communicator);
        }
        {
          Table<2, DoFTools::Coupling> coupling(2 * spacedim, 2 * spacedim);
          for (unsigned int c = 0; c < spacedim; ++c)
            for (unsigned int d = 0; d < spacedim; ++d)
              coupling[c][d] = DoFTools::none;
          for (unsigned int c = spacedim; c < 2 * spacedim; ++c)
            for (unsigned int d = spacedim; d < 2 * spacedim; ++d)
              coupling[c][d] = (c == d) ? DoFTools::always : DoFTools::none;

          BlockDynamicSparsityPattern dsp(problem.solid_dofs_per_block,
                                          problem.solid_dofs_per_block);
          DoFTools::make_sparsity_pattern(
            problem.solid_dh, coupling, dsp, problem.solid_constraints, false);
          SparsityTools::distribute_sparsity_pattern(
            dsp,
            locally_owned_solid_dofs_per_processor,
            problem.mpi_communicator,
            locally_relevant_solid_dofs);
          problem.solid_preconditioner.reinit(problem.solid_owned_dofs,
                                              dsp,
                                              problem.mpi_communicator);
        }
      }

      // The previous displacement w_old is zero (setup_dofs() reinitializes it
      // to zero), which is exactly what the manufactured multiplier forcing
      // assumes.
      const double alpha = 1.0;

      // call assemble functions from the problem
      problem.assemble_elasticity_system(alpha);
      problem.assemble_elasticity_rhs(alpha);
      solve_elasticity(problem);

      // L2 error of the displacement components only (no multiplier).
      const ComponentSelectFunction<spacedim> displacement_mask(
        std::pair<unsigned int, unsigned int>(0, spacedim), 2 * spacedim);

      Vector<double> cellwise_error(
        problem.solid_dh.get_triangulation().n_active_cells());
      VectorTools::integrate_difference(problem.solid_dh,
                                        problem.solid_locally_relevant_solution,
                                        exact_solution,
                                        cellwise_error,
                                        QGauss<spacedim>(
                                          problem.solid_fe->degree + 2),
                                        VectorTools::L2_norm,
                                        &displacement_mask);

      const double l2_error =
        VectorTools::compute_global_error(problem.solid_dh.get_triangulation(),
                                          cellwise_error,
                                          VectorTools::L2_norm);

      const double h = 2.0 / std::pow(2.0, level);
      h_values.push_back(h);
      l2_errors.push_back(l2_error);

      convergence_table.add_value("level", level);
      convergence_table.add_value(
        "cells", problem.solid_dh.get_triangulation().n_active_cells());
      convergence_table.add_value("dofs", problem.solid_dh.n_dofs());
      convergence_table.add_value("h", h);
      convergence_table.add_value("L2_displacement", l2_error);
    }

  convergence_table.set_precision("h", 4);
  convergence_table.set_scientific("h", true);
  convergence_table.set_precision("L2_displacement", 6);
  convergence_table.set_scientific("L2_displacement", true);
  convergence_table.evaluate_convergence_rates(
    "L2_displacement", ConvergenceTable::reduction_rate_log2);

  std::cout << "\nElasticity manufactured-solution convergence:\n";
  convergence_table.write_text(std::cout);
  std::cout << std::endl;

  // We do the same we did for the Stokes test. Now with Q_1 displacement, we
  // expect an L2 convergence rate of 2. Require the observed rate on each
  // successive mesh pair to exceed displacement_degree + 0.5.
  for (unsigned int i = 1; i < levels.size(); ++i)
    {
      const double rate = std::log(l2_errors[i - 1] / l2_errors[i]) /
                          std::log(h_values[i - 1] / h_values[i]);
      EXPECT_GT(rate, 1.5)
        << "L2 displacement convergence rate between refinement levels "
        << levels[i - 1] << " and " << levels[i] << " is only " << rate
        << " (errors " << l2_errors[i - 1] << " -> " << l2_errors[i] << ").";
    }
}



namespace Step80
{
  // Pointwise internal force of the compressible neo-Hookean law used in the
  // test below to check the tangent formulas by finite differences:
  //   N(F) = (P(F), G),  P(F) = mu (F - F^{-T}) + lambda ln(J) F^{-T}.
  template <int dim>
  double pointwise_neohookean_force(const Tensor<2, dim> &F,
                                    const Tensor<2, dim> &G,
                                    const double          mu,
                                    const double          lambda)
  {
    using namespace dealii;

    const Tensor<2, dim> Q   = transpose(invert(F));
    const double         lnJ = std::log(determinant(F));
    return mu * scalar_product(F, G) +
           (lambda * lnJ - mu) * scalar_product(Q, G);
  }



  // Pointwise internal force of the exponential law
  // P(F) = mu exp(lambda (F:F - dim)) F.
  template <int dim>
  double pointwise_exponential_force(const Tensor<2, dim> &F,
                                     const Tensor<2, dim> &G,
                                     const double          mu,
                                     const double          lambda)
  {
    using namespace dealii;

    const double g =
      mu * std::exp(lambda * (scalar_product(F, F) - static_cast<double>(dim)));
    return g * scalar_product(F, G);
  }



  // Pointwise verification of the stress and material tangent of the two
  // hyperelastic laws implemented in the tutorial, against finite differences
  // of the internal force. The test is run both in 2D and in 3D.
  template <int dim>
  void check_hyperelastic_tangent_pointwise()
  {
    using namespace dealii;

    const double mu      = 1.3;
    const double lambda  = 2.7;
    const double epsilon = 1e-6;

    // Gradients used to contract the tangent. Keep H small so that the
    // configuration F = I + H used by the right-hand side entry stays
    // invertible.
    Tensor<2, dim> G_a, G_b, H;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        {
          G_a[i][j] = 0.3 + 0.2 * i - 0.1 * j;
          G_b[i][j] = -0.2 + 0.15 * j + 0.25 * (i == j);
          H[i][j]   = 0.03 * (i + 1) - 0.02 * j;
        }

    // Data evaluated from the displacement gradient H: F_data = I + H.
    const HyperelasticPointData<dim> data_nh =
      evaluate_hyperelastic_point_data(H, true, "neo_hookean", mu, lambda);
    Tensor<2, dim> F_nh;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        F_nh[i][j] = (i == j ? 1. : 0.) + H[i][j];

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        EXPECT_NEAR(data_nh.F[i][j], F_nh[i][j], 1e-14);
    EXPECT_NEAR(data_nh.ln_det_F, std::log(determinant(F_nh)), 1e-14);

    // Tangent entry against central finite differences of the internal force:
    // (dP/dF[G_b], G_a) ~ (N(F + eps G_b) - N(F - eps G_b)) / (2 eps).
    const Tensor<2, dim> F_plus  = F_nh + epsilon * G_b;
    const Tensor<2, dim> F_minus = F_nh - epsilon * G_b;
    const double         fd_tangent =
      (pointwise_neohookean_force(F_plus, G_a, mu, lambda) -
       pointwise_neohookean_force(F_minus, G_a, mu, lambda)) /
      (2. * epsilon);

    const double tangent =
      neohookean_tangent_entry(G_a, G_b, data_nh, mu, lambda);
    EXPECT_NEAR(tangent, fd_tangent, 1e-7);

    // Right-hand side entry: (K_s(w) w)_i - N_i, where the tangent is applied
    // to H = grad(w) and the internal force is evaluated at F = I + H.
    const double Ks_w = neohookean_tangent_entry(G_a, H, data_nh, mu, lambda);
    const double N_i  = pointwise_neohookean_force(F_nh, G_a, mu, lambda);
    EXPECT_NEAR(neohookean_rhs_entry(G_a, data_nh, H, mu, lambda),
                Ks_w - N_i,
                1e-14);

    // At F = I the neo-Hookean tangent must reduce to the linear elastic
    // operator 2 mu sym(G_a) : sym(G_b) + lambda div(G_a) div(G_b), and the
    // internal force must vanish.
    {
      const Tensor<2, dim>             H_zero;
      const HyperelasticPointData<dim> data_identity =
        evaluate_hyperelastic_point_data(
          H_zero, true, "neo_hookean", mu, lambda);
      const double linear_entry =
        2. * mu * scalar_product(symmetrize(G_a), symmetrize(G_b)) +
        lambda * trace(G_a) * trace(G_b);
      EXPECT_NEAR(neohookean_tangent_entry(G_a, G_b, data_identity, mu, lambda),
                  linear_entry,
                  1e-12);
      EXPECT_NEAR(pointwise_neohookean_force(data_identity.F, G_a, mu, lambda),
                  0.,
                  1e-12);
    }

    // The same checks for the exponential law.
    const HyperelasticPointData<dim> data_exp =
      evaluate_hyperelastic_point_data(H, true, "exponential", mu, lambda);
    const Tensor<2, dim> F_exp = F_nh;
    const double         fd_tangent_exp =
      (pointwise_exponential_force(F_exp + epsilon * G_b, G_a, mu, lambda) -
       pointwise_exponential_force(F_exp - epsilon * G_b, G_a, mu, lambda)) /
      (2. * epsilon);
    EXPECT_NEAR(exponential_tangent_entry(G_a, G_b, data_exp, lambda),
                fd_tangent_exp,
                1e-7);
    const double Ks_w_exp = exponential_tangent_entry(G_a, H, data_exp, lambda);
    const double N_i_exp  = pointwise_exponential_force(F_exp, G_a, mu, lambda);
    EXPECT_NEAR(exponential_rhs_entry(G_a, data_exp, H, lambda),
                Ks_w_exp - N_i_exp,
                1e-14);
  }



} // namespace Step80



namespace Step80
{
  // Independent assembly of the nonlinear internal force
  //   N_i(w) = int_B P(F) : grad(phi_i) dX
  // on the solid mesh, with F = I + grad(w) and P(F) the first
  // Piola-Kirchhoff stress of the law selected in par. The stress formulas
  // are written out directly here, i.e. independently of the entry functions
  // used by the tutorial assembly, so that the assembled tangent stiffness
  // and right-hand side can be checked against a second implementation.
  template <int dim, int spacedim>
  LA::MPI::Vector assemble_internal_force(
    const NavierStokesImmersedProblem<dim, spacedim> &problem)
  {
    using namespace dealii;

    const auto &par = problem.par;

    FEValues<spacedim> fe_values(*problem.solid_fe,
                                 problem.solid_quadrature,
                                 update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = problem.solid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = problem.solid_quadrature.size();

    std::vector<Tensor<2, spacedim>>     grad_w(n_q_points);
    Vector<double>                       cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    LA::MPI::BlockVector force(problem.solid_solution);
    force = 0;

    for (const auto &cell : problem.solid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[problem.displacement].get_function_gradients(
            problem.solid_locally_relevant_solution, grad_w);

          cell_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Deformation gradient F = I + grad(w).
              Tensor<2, spacedim> F = grad_w[q];
              for (unsigned int d = 0; d < spacedim; ++d)
                F[d][d] += 1.;

              // First Piola-Kirchhoff stress of the selected law.
              Tensor<2, spacedim> P;
              if (par.nonlinear_solid_law == "neo_hookean")
                {
                  const Tensor<2, spacedim> Q   = transpose(invert(F));
                  const double              lnJ = std::log(determinant(F));
                  P = par.lame_mu * (F - Q) + par.lame_lambda * lnJ * Q;
                }
              else
                {
                  const double g =
                    par.lame_mu *
                    std::exp(par.lame_lambda * (scalar_product(F, F) -
                                                static_cast<double>(spacedim)));
                  P = g * F;
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) +=
                  scalar_product(P,
                                 fe_values[problem.displacement].gradient(i,
                                                                          q)) *
                  fe_values.JxW(q);
            }

          cell->get_dof_indices(local_dof_indices);
          problem.solid_constraints.distribute_local_to_global(
            cell_rhs, local_dof_indices, force);
        }

    force.compress(VectorOperation::add);
    return LA::MPI::Vector(force.block(0));
  }
} // namespace Step80



TEST(Step80, HyperelasticStressAndTangent)
{
  using namespace Step80;

  check_hyperelastic_tangent_pointwise<2>();
  check_hyperelastic_tangent_pointwise<3>();
}



// End-to-end check of the assembled nonlinear solid terms for both
// hyperelastic laws, against an independent assembly of the internal force:
//
//  1. The assembled tangent stiffness K_s(w) applied to a test direction x
//     is compared with the centered finite difference of the independent
//     internal force, (N(w + eps x) - N(w - eps x)) / (2 eps).
//  2. With vanishing previous displacement and body force, the solid
//     right-hand side assembled by assemble_elasticity_rhs() is
//     R(w) = K_s(w) w - N(w), so it must match the matrix-vector product
//     K_s(w) w minus the independently assembled N(w).
TEST(Step80, NonlinearSolidAssemblyTangentMatchesRightHandSide)
{
  using namespace dealii;
  using namespace Step80;

  ParameterAcceptor::clear();

  constexpr int dim      = 2;
  constexpr int spacedim = 2;

  std::ostringstream prm;
  prm << "subsection Navier-Stokes Immersed Problem\n"
      << "  subsection Grid generation\n"
      << "    set Fluid grid generator           = hyper_cube\n"
      << "    set Fluid grid generator arguments = -1: 1: false\n"
      << "    set Initial fluid refinement       = 1\n"
      << "    set Solid grid generator           = hyper_cube\n"
      << "    set Solid grid generator arguments = -1: 1: false\n"
      << "    set Initial solid refinement       = 3\n"
      << "  end\n"
      << "  subsection Physical properties\n"
      << "    set Lame mu               = 2.5\n"
      << "    set Lame lambda           = 3.7\n"
      << "  end\n"
      << "end\n";

  const std::string prm_filename =
    (std::filesystem::temp_directory_path() / "step-80-nonlinear-tangent.prm")
      .string();
  {
    std::ofstream out(prm_filename);
    out << prm.str();
  }

  const auto [file_dim, file_spacedim] =
    get_dimension_and_spacedimension(prm_filename);
  ASSERT_EQ(file_dim, static_cast<unsigned int>(dim));
  ASSERT_EQ(file_spacedim, static_cast<unsigned int>(spacedim));

  NavierStokesImmersedProblemParameters<dim, spacedim> par;
  ParameterAcceptor::initialize(prm_filename);

  const double alpha = 1.0;
  const double eps   = 1e-6;

  for (const std::string law : {"neo_hookean", "exponential"})
    {
      par.use_nonlinear_solid_model = true;
      par.nonlinear_solid_law       = law;

      NavierStokesImmersedProblem<dim, spacedim> problem(par);
      problem.make_grid();
      problem.initial_setup();
      problem.setup_dofs();

      // Deterministic test displacement w on the solid displacement block
      // (the multiplier block stays at zero), small enough that the
      // deformation gradient F = I + grad(w) remains invertible, and a test
      // direction x for the finite difference.
      LA::MPI::Vector w(problem.solid_owned_dofs[0], problem.mpi_communicator);
      LA::MPI::Vector x(problem.solid_owned_dofs[0], problem.mpi_communicator);
      for (auto i : w.locally_owned_elements())
        {
          w[i] = 0.02 * std::sin(0.37 * static_cast<double>(i) + 0.5);
          x[i] = 0.01 * std::cos(0.23 * static_cast<double>(i) + 0.1);
        }

      auto set_displacement = [&](const LA::MPI::BlockVector &v) {
        problem.solid_solution                  = v;
        problem.solid_locally_relevant_solution = v;
      };

      // 1) Tangent stiffness at w against a finite difference of the
      // independent internal force.
      LA::MPI::BlockVector w_plus(problem.solid_solution);
      LA::MPI::BlockVector w_minus(problem.solid_solution);
      w_plus           = 0;
      w_minus          = 0;
      w_plus.block(0)  = w;
      w_minus.block(0) = w;
      w_plus.block(0).add(eps, x);
      w_minus.block(0).add(-eps, x);

      set_displacement(w_plus);
      const LA::MPI::Vector N_plus = assemble_internal_force(problem);
      set_displacement(w_minus);
      const LA::MPI::Vector N_minus = assemble_internal_force(problem);

      LA::MPI::Vector fd(problem.solid_owned_dofs[0], problem.mpi_communicator);
      fd = N_plus;
      fd.add(-1., N_minus);
      fd *= 1. / (2. * eps);

      LA::MPI::BlockVector w_0(problem.solid_solution);
      w_0          = 0;
      w_0.block(0) = w;
      set_displacement(w_0);
      problem.assemble_elasticity_system(alpha);
      const auto &K = problem.solid_matrix.block(0, 0);

      LA::MPI::Vector Kx(problem.solid_owned_dofs[0], problem.mpi_communicator);
      K.vmult(Kx, x);
      Kx.add(-1., fd);

      const double tangent_error = Kx.l2_norm() / x.l2_norm();
      EXPECT_LT(tangent_error, 1e-6)
        << "Assembled tangent stiffness and finite difference of the "
        << "internal force disagree for law '" << law << "': relative error "
        << tangent_error;

      // 2) Assembled right-hand side R(w) = K_s(w) w - N(w) against the
      // matrix-vector product and the independent internal force.
      problem.assemble_elasticity_rhs(alpha);
      const LA::MPI::Vector R = problem.solid_system_rhs.block(0);
      const LA::MPI::Vector N = assemble_internal_force(problem);

      LA::MPI::Vector Kw(problem.solid_owned_dofs[0], problem.mpi_communicator);
      K.vmult(Kw, w);

      LA::MPI::Vector residual(problem.solid_owned_dofs[0],
                               problem.mpi_communicator);
      residual = R;
      residual.add(-1., Kw);
      residual.add(1., N);

      const double rhs_error = residual.l2_norm() / w.l2_norm();
      EXPECT_LT(rhs_error, 1e-9)
        << "Assembled right-hand side does not match K_s(w) w - N(w) for law '"
        << law << "': relative error " << rhs_error;
    }
}

int main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
