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

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "step-80.cc"

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

    std::array<std::array<LinOp, 2>, 2> system_array = {
      {{{A11_aug, Bt}},  // vel
       {{B, Z6}}}};      // pres

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

    SolverControl solver_control(par.outer_max_iterations,
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
} // namespace Step80

TEST(Step80, GetDimensionAndSpacedimensionReadsParameterHandler) {
  dealii::ParameterHandler prm;

  prm.declare_entry("dimension", "2", dealii::Patterns::Integer(2, 3));
  prm.declare_entry("space dimension", "2", dealii::Patterns::Integer(2, 3));
  prm.set("dimension", "3");
  prm.set("space dimension", "3");

  const auto [dim, spacedim] = Step80::get_dimension_and_spacedimension(prm);

  EXPECT_EQ(dim, 3U);
  EXPECT_EQ(spacedim, 3U);
}

TEST(Step80, MPIIsInitialized) {
  int initialized = 0;
  ASSERT_EQ(MPI_Initialized(&initialized), MPI_SUCCESS);
  EXPECT_NE(initialized, 0);

  EXPECT_EQ(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 1U);
}

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
TEST(Step80, NavierStokesManufacturedSolutionConvergence) {
  using namespace dealii;
  using namespace Step80;

  constexpr int dim = 2;
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
  std::vector<double> h_values;
  std::vector<double> l2_errors;

  ConvergenceTable convergence_table;

  for (const unsigned int level : levels) {
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
    for (unsigned int it = 0; it < 50 && change > 1e-9; ++it) {
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
    VectorTools::integrate_difference(
        problem.fluid_dh, problem.fluid_locally_relevant_solution,
        par.navier_stokes_analytical_solution, cellwise_error,
        QGauss<spacedim>(problem.fluid_fe->degree + 2), VectorTools::L2_norm,
        &velocity_mask);

    const double l2_error =
        VectorTools::compute_global_error(problem.fluid_dh.get_triangulation(),
                                          cellwise_error, VectorTools::L2_norm);

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
  for (unsigned int i = 1; i < levels.size(); ++i) {
    const double rate = std::log(l2_errors[i - 1] / l2_errors[i]) /
                        std::log(h_values[i - 1] / h_values[i]);
    EXPECT_GT(rate, 2.5)
        << "L2 velocity convergence rate between refinement levels "
        << levels[i - 1] << " and " << levels[i] << " is only " << rate
        << " (errors " << l2_errors[i - 1] << " -> " << l2_errors[i] << ").";
  }
}

int main(int argc, char *argv[]) {
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
