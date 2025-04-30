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



// Check LinearAlgebra::TpetraWrappers::NOXSolver by solving f(x) = atan(x)-0.5
// = 0, starting at x=10. This corresponds also to what the sundials/kinsol_06
// test checks.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/trilinos/nox.h>

// as reference solution
#include <deal.II/lac/trilinos_tpetra_precondition.h>
#include <deal.II/lac/trilinos_tpetra_solver.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  using Number     = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  // set up solver control
  const unsigned int n_max_iterations  = 100;
  const double       abs_tolerance     = 1e-9;
  const double       rel_tolerance     = 1e-5;
  const double       lin_rel_tolerance = 1e-3;

  LinearAlgebra::TpetraWrappers::NOXSolver<VectorType>::AdditionalData
    additional_data(n_max_iterations, abs_tolerance, rel_tolerance);

  // set up parameters
  Teuchos::RCP<Teuchos::ParameterList> non_linear_parameters =
    Teuchos::rcp(new Teuchos::ParameterList);

  non_linear_parameters->set("Nonlinear Solver", "Line Search Based");
  non_linear_parameters->sublist("Printing").set("Output Information", 15);
  non_linear_parameters->sublist("Direction").set("Method", "Newton");
  non_linear_parameters->sublist("Direction")
    .sublist("Newton")
    .sublist("Linear Solver")
    .set("Tolerance", lin_rel_tolerance);
  non_linear_parameters->sublist("Line Search").set("Method", "Polynomial");


  // set up solver
  LinearAlgebra::TpetraWrappers::NOXSolver<VectorType> solver(
    additional_data, non_linear_parameters);

  // ... helper functions
  double J = 0.0;

  solver.residual = [](const auto &src, auto &dst) {
    // compute residual
    deallog << "Evaluating residual at u=" << src[0] << std::endl;
    dst[0] = std::atan(src[0]) - 0.5;
  };

  solver.setup_jacobian = [&](const auto &src) {
    // compute Jacobian
    deallog << "Setting up Jacobian at u=" << src[0] << std::endl;
    J = 1. / (1 + src[0] * src[0]);
  };

  solver.apply_jacobian = [&](const auto &src, auto &dst) {
    // solve with Jacobian
    dst[0] = src[0] * J;
  };

  solver.solve_with_jacobian = [&](const auto &src, auto &dst, const auto) {
    // solve with Jacobian
    dst[0] = src[0] / J;
  };

  // initial guess
  VectorType solution(1);
  solution[0] = 10.0;

  // solve with the given initial guess
  solver.solve(solution);

  deallog << "The solution is: " << solution[0] << std::endl;
}
