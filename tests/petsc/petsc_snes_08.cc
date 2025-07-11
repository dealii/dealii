// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/utilities.h>

#include <deal.II/lac/petsc_vector.h>

#include <deal.II/numerics/nonlinear.h>

#include "../tests.h"


// Like the sundials/kinsol_06_v3 test, but for SNES. After a few
// tries with the given residual/Jacobian, SNES gives up and
// terminates because it can't find a solution. We need to make sure
// that we throw a catchable exception.
//
// This testcase is a variation of the kinsol_06 test, modified by
// Simon Wiesheier, and posted on the mailing list. Then further
// adapted by Wolfgang Bangerth.


int
main(int argc, char *argv[])
{
  initlog();

  using VectorType = PETScWrappers::MPI::Vector;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  // Size of the problem
  const unsigned int N = 1;

  NonlinearSolverSelector<VectorType>::AdditionalData additional_data;
  additional_data.solver_type =
    NonlinearSolverSelector<VectorType>::AdditionalData::SolverType::petsc_snes;
  additional_data.maximum_non_linear_iterations = 1;

  NonlinearSolverSelector<VectorType> nonlinear_solver(additional_data);

  nonlinear_solver.reinit_vector = [N](VectorType &v) {
    v.reinit(MPI_COMM_WORLD, N, N);
  };

  nonlinear_solver.residual = [&](const VectorType &u, VectorType &F) {
    // Count how many times this function has been called.
    static int count = 0;
    deallog << "Computing residual for the " << count + 1
            << "th time, at u=" << u[0] << std::endl;
    if ((u[0] < -10) || (u[0] > 20) /* || (count > 0) */)
      {
        deallog << "Reporting recoverable failure." << std::endl;
        throw RecoverableUserCallbackError();
      }
    count++;


    F.reinit(u);
    F[0] = std::atan(u[0]) - 0.5;
    F.compress(VectorOperation::values::insert);
  };

  double J_inverse;

  nonlinear_solver.setup_jacobian = [&J_inverse](const VectorType &u) {
    deallog << "Setting up Jacobian system at u=" << u[0] << std::endl;

    const double J = 1. / (1 + u[0] * u[0]);
    J_inverse      = 1. / J;
  };


  nonlinear_solver.solve_with_jacobian =
    [&](const VectorType &rhs, VectorType &dst, double) {
      dst[0] = J_inverse * rhs[0];
      dst.compress(VectorOperation::values::insert);
    };


  VectorType v(MPI_COMM_WORLD, N, N);
  v[0] = 10;
  v.compress(VectorOperation::values::insert);

  try
    {
      nonlinear_solver.solve(v);
    }
  catch (const ExceptionBase &e)
    {
      deallog
        << "Nonlinear Solver threw an exception with the following message:"
        << std::endl;
      e.print_info(deallog.get_file_stream());
    }
}
