//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2021 by the deal.II authors
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

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"

// Solve a nonlinear system but provide only residual function. KINSOL
// then uses its internal solvers which are based on a
// finite-difference approximation to the Jacobian and a direct
// solver.
//
// Compared to the _01 test, this is simply a more complicated function:
// We solve the nonlinear problem
//
//   F(u) = 0
//
// with a 2-dimensional vector u and where
//
//   F(u) = [  cos(u1 + u2) - 1   ]              -> u1=-u2
//          [  sin(u1 - u2)       ]              -> u1=u2
//
// In other words, we need to find the solution u1=u2=0.

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  using VectorType = Vector<double>;

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/kinsol_01.prm");
  prm.parse_input(ifile);

  // Size of the problem
  unsigned int N = 2;

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [](const VectorType &u, VectorType &F) -> int {
    F(0) = std::cos(u[0] + u[1]) - 1;
    F(1) = std::sin(u[0] - u[1]);
    return 0;
  };


  kinsol.iteration_function = [](const VectorType &u, VectorType &F) -> int {
    // We want a Newton-type scheme, not a fixed point iteration. So we
    // shouldn't get into this function.
    std::abort();

    // But if anyone wanted to see how it would look like:
    F(0) = std::cos(u[0] + u[1]) - 1 - u[0];
    F(1) = std::sin(u[0] - u[1]) - u[1];
    return 0;
  };

  VectorType v(N);
  v(0) = 0.5;
  v(1) = 1.234;

  auto niter = kinsol.solve(v);
  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
