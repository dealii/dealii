// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Check LinearAlgebra::TpetraWrappers::NOXSolver. This test solves the same
// problem as nox_solver_03, but it throws an exception in solve_with_jacobian.
// This callback is special, since with default settings NOX tries to perform a
// recovery step once an exception has been raised there. We need to treat such
// situations accordingly. This test verifies that the recovery functionality of
// NOX works.


#include <deal.II/base/mpi.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/trilinos/nox.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  using Number     = double;
  using VectorType = Vector<Number>;

  // set up solver control
  const unsigned int n_max_iterations = 10;
  const double       abs_tolerance    = 1e-9;
  const double       rel_tolerance    = 1e-5;

  LinearAlgebra::TpetraWrappers::NOXSolver<VectorType>::AdditionalData
    additional_data(n_max_iterations, abs_tolerance, rel_tolerance);

  // set up parameters
  Teuchos::RCP<Teuchos::ParameterList> non_linear_parameters =
    Teuchos::rcp(new Teuchos::ParameterList);

  non_linear_parameters->set("Nonlinear Solver", "Line Search Based");
  non_linear_parameters->sublist("Printing").set("Output Information", 15);
  non_linear_parameters->sublist("Direction").set("Method", "Newton");
  non_linear_parameters->sublist("Line Search").set("Method", "Polynomial");


  // set up solver
  LinearAlgebra::TpetraWrappers::NOXSolver<VectorType> solver(
    additional_data, non_linear_parameters);

  // ... helper functions
  solver.residual = [](const VectorType &u, VectorType &F) {
    deallog << "Evaluating the solution at u=(" << u[0] << ',' << u[1] << ')'
            << std::endl;

    F(0) = std::cos(u[0] + u[1]) - 1 + 2 * u[0];
    F(1) = std::sin(u[0] - u[1]) + 2 * u[1];
  };

  FullMatrix<double> J(2, 2);
  FullMatrix<double> J_inverse(2, 2);

  solver.setup_jacobian = [&J, &J_inverse](const VectorType &u) {
    // We don't do any kind of set-up in this program, but we can at least
    // say that we're here
    deallog << "Setting up Jacobian system at u=(" << u[0] << ',' << u[1] << ')'
            << std::endl;

    J(0, 0) = -std::sin(u[0] + u[1]) + 2;
    J(0, 1) = -std::sin(u[0] + u[1]);
    J(1, 0) = std::cos(u[0] - u[1]);
    J(1, 1) = -std::cos(u[0] - u[1]) + 2;

    J_inverse.invert(J);
  };

  solver.apply_jacobian = [&](const VectorType &src, VectorType &dst) {
    J.vmult(dst, src);
  };

  unsigned int counter = 0;

  solver.solve_with_jacobian = [&J_inverse,
                                &counter](const VectorType &rhs,
                                          VectorType       &dst,
                                          const double /*tolerance*/) {
    deallog << "Solving Jacobian system with rhs=(" << rhs[0] << ',' << rhs[1]
            << ')' << std::endl;

    J_inverse.vmult(dst, rhs);

    counter++;

    // Fail linear solver step but NOX should perform a recovery step
    if (counter == 2)
      throw RecoverableUserCallbackError(
        "Recoverable failure in solve_with_jacobian.");

    // Throw an irrecoverable exception then
    if (counter == 5)
      throw ExcMessage("Irrecoverable failure in solve_with_jacobian.");
  };

  // initial guess
  const unsigned int N = 2;
  VectorType         solution(N);
  solution[0] = 0.5;
  solution[1] = 1.234;

  // solve with the given initial guess
  try
    {
      auto niter = solver.solve(solution);
      solution.print(deallog.get_file_stream());
      deallog << "Converged in " << niter << " iterations." << std::endl;
    }
  catch (const std::exception &exc)
    {
      deallog << "Caught an irrecoverable exception in a callback:" << std::endl
              << exc.what() << std::endl;
    }
  catch (const char *s)
    {
      deallog << "Trilinos threw its own exception of type char*:" << std::endl
              << s << std::endl;
    }
}
