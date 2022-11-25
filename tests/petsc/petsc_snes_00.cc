//-----------------------------------------------------------
//
//    Copyright (C) 2022 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------
//
// Author: Stefano Zampini, King Abdullah University of Science and Technology.

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_snes.h>
#include <deal.II/lac/petsc_vector.h>

#include <cmath>

#include "../tests.h"

/**
 * Solves the nonlinear system of equations
 *
 * (x - y^3 + 1)^3 - y^3 = 0
 * x + 2y - 3 = 0
 *
 * using the PETScWrappers::NonlinearSolver class
 * that interfaces PETSc SNES solver object.
 *
 * Here we test a pure SNES interface and a deal.ii approach
 * when users take control of solving the linear systems.
 */
using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MatrixBase;
using Solver     = PETScWrappers::NonlinearSolver<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream out("output");

  PETScWrappers::NonlinearSolverData data;
  ParameterHandler                   prm;

  data.add_parameters(prm);
  out << "# Default Parameters" << std::endl;
  prm.print_parameters(out, ParameterHandler::ShortText);

  std::ifstream ifile(SOURCE_DIR "/petsc_snes_00_in.prm");
  prm.parse_input(ifile);

  out << "# Testing Parameters" << std::endl;
  prm.print_parameters(out, ParameterHandler::ShortText);

  {
    out << "# Test PETSc interface (FD)" << std::endl;

    // Create and customize the solver
    Solver solver(data);

    // Here we pass the callback for the function evaluation
    // If the solver needs the Jacobian, PETSc will compute it with
    // expensive finite-differencing of the user function.
    solver.residual = [&](const VectorType &X, VectorType &F) -> int {
      auto x = X[0];
      auto y = X[1];
      F(0)   = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F(1)   = x + 2 * y - 3;
      F.compress(VectorOperation::insert);
      return 0;
    };

    // Create and initialize solution vector
    VectorType x(MPI_COMM_SELF, 2, 2);
    x = 0.0;

    // Solve the nonlinear system (defaults to Newton with cubic-backtracking)
    auto nit = solver.solve(x);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }

  {
    out << "# Test PETSc interface" << std::endl;

    Solver solver(data);

    solver.residual = [&](const VectorType &X, VectorType &F) -> int {
      auto x = X[0];
      auto y = X[1];
      F(0)   = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F(1)   = x + 2 * y - 3;
      F.compress(VectorOperation::insert);
      return 0;
    };

    // Here we use the Jacobian callback following the PETSc style,
    // where matrices are used as output arguments.
    // This test uses A == P, so we do not do anything on A
    // The intermediate callbacks used to pass from PETSc to deal.ii
    // will handle the case of users requesting Jacobian-free Newton
    // Krylov (i.e. using -snes_mf_operator)
    solver.jacobian =
      [&](const VectorType &X, MatrixType &A, MatrixType &P) -> int {
      auto x    = X[0];
      auto y    = X[1];
      auto f0_x = 3 * std::pow(x - std::pow(y, 3) + 1, 2);
      auto f0_y = -9 * std::pow(x - std::pow(y, 3) + 1, 2) * std::pow(y, 2) -
                  3 * std::pow(y, 2);
      P.set(0, 0, f0_x);
      P.set(0, 1, f0_y);
      P.set(1, 0, 1);
      P.set(1, 1, 2);
      P.compress(VectorOperation::insert);
      return 0;
    };

    VectorType x(MPI_COMM_SELF, 2, 2);
    x = 0.0;

    // Solve the nonlinear system without specifying matrices
    // This will use the default matrix type in PETSc for SNES
    // (MATDENSE), which is ok for our test.
    // See petsc_snes_01.cc for a case where we pass the matrix
    auto nit = solver.solve(x);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }

  {
    out << "# Test user interface" << std::endl;

    Solver solver(data);

    solver.residual = [&](const VectorType &X, VectorType &F) -> int {
      auto x = X[0];
      auto y = X[1];
      F(0)   = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F(1)   = x + 2 * y - 3;
      F.compress(VectorOperation::insert);
      return 0;
    };

    // When users want to be in full control of the linear system
    // solves, they need to use setup_jacobian and solve_for_jacobian_system
    // For example, in this case we compute the inverse of the Jacobian
    // during setup_jacobian and we use it in the solve phase
    FullMatrix<double> Jinv(2, 2);

    solver.setup_jacobian = [&](const VectorType &X) -> int {
      auto x    = X[0];
      auto y    = X[1];
      auto f0_x = 3 * std::pow(x - std::pow(y, 3) + 1, 2);
      auto f0_y = -9 * std::pow(x - std::pow(y, 3) + 1, 2) * std::pow(y, 2) -
                  3 * std::pow(y, 2);
      FullMatrix<double> J(2, 2);
      J(0, 0) = f0_x;
      J(0, 1) = f0_y;
      J(1, 0) = 1;
      J(1, 1) = 2;
      Jinv.invert(J);
      return 0;
    };

    // Solve phase. By default, PETSc will use this callback as a preconditioner
    // within a preconditioner only Krylov solve. Other Krylov
    // solvers can still be used in a Jacobian-free way and selected at command
    // line or within user code.
    solver.solve_for_jacobian_system = [&](const VectorType &src,
                                           VectorType &      dst) -> int {
      dst(0) = Jinv(0, 0) * src(0) + Jinv(0, 1) * src(1);
      dst(1) = Jinv(1, 0) * src(0) + Jinv(1, 1) * src(1);
      dst.compress(VectorOperation::insert);
      return 0;
    };

    VectorType x(MPI_COMM_SELF, 2, 2);
    x = 0.0;

    // We don't need to pass matrices, since in this case we don't need them
    auto nit = solver.solve(x);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }

  {
    out << "# Test user interface (affine)" << std::endl;

    Solver solver(data);

    solver.residual = [&](const VectorType &X, VectorType &F) -> int {
      auto x = X[0];
      auto y = X[1];
      F(0)   = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F(1)   = x + 2 * y;
      F.compress(VectorOperation::insert);
      return 0;
    };

    FullMatrix<double> Jinv(2, 2);

    solver.setup_jacobian = [&](const VectorType &X) -> int {
      auto x    = X[0];
      auto y    = X[1];
      auto f0_x = 3 * std::pow(x - std::pow(y, 3) + 1, 2);
      auto f0_y = -9 * std::pow(x - std::pow(y, 3) + 1, 2) * std::pow(y, 2) -
                  3 * std::pow(y, 2);
      FullMatrix<double> J(2, 2);
      J(0, 0) = f0_x;
      J(0, 1) = f0_y;
      J(1, 0) = 1;
      J(1, 1) = 2;
      Jinv.invert(J);
      return 0;
    };

    solver.solve_for_jacobian_system = [&](const VectorType &src,
                                           VectorType &      dst) -> int {
      dst(0) = Jinv(0, 0) * src(0) + Jinv(0, 1) * src(1);
      dst(1) = Jinv(1, 0) * src(0) + Jinv(1, 1) * src(1);
      dst.compress(VectorOperation::insert);
      return 0;
    };

    VectorType x(MPI_COMM_SELF, 2, 2);
    x = 0.0;

    // Here we test the affine interface for solving F(u) = b
    // Create the b = [0,3] vector
    VectorType rhs(x);
    rhs[1] = 3.0;
    rhs.compress(VectorOperation::insert);

    auto nit = solver.solve(x, rhs);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }
  return 0;
}
