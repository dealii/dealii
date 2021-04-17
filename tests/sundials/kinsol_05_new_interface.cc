//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
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

// Solve a nonlinear system.
//
// Similar to the _02 test, but we're now actually providing a solver with
// the Jacobian matrix. For the current case,
//
//   F(u) = [  cos(u1 + u2) - 1 + 2*u1  ]
//          [  sin(u1 - u2)     + 2*u2   ]
//
// the Jacobian is the 2x2 matrix
//
//   J(u) = [ -sin(u1 + u2) + 2      -sin(u1 + u2)]
//          [  cos(u1 - u2)          -cos(u1 - u2) + 2]
//
// The addition of the +2u_i to the function F does not move the solution
// (it is still u=0) but it makes sure that the Jacobian at the solution
// remains non-singular
//
// The _03 test case has a flaw in that starting at SUNDIALS 4.x, the solve
// function no longer receives the current 'u' vector. This means that
// one can't compute a proper Jacobian in that function because we don't
// know what to linearize around. This test fixes this by computing
// the Jacobian in the setup function.
//
// The _04 test fixes this by setting up the Jacobian in the setup
// function, which has the current iterate. But it converges slowly
// because it doesn't update the Jacobian very often. This test
// finally updates it in every iteration.

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

  // Update the Jacobian in each iteration:
  data.maximum_setup_calls = 1;


  // Size of the problem
  unsigned int N = 2;

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [](const VectorType &u, VectorType &F) -> int {
    deallog << "Evaluating the solution at u=(" << u[0] << ',' << u[1] << ")"
            << std::endl;

    F(0) = std::cos(u[0] + u[1]) - 1 + 2 * u[0];
    F(1) = std::sin(u[0] - u[1]) + 2 * u[1];
    return 0;
  };


  kinsol.iteration_function = [](const VectorType &u, VectorType &F) -> int {
    // We want a Newton-type scheme, not a fixed point iteration. So we
    // shouldn't get into this function.
    std::abort();

    // But if anyone wanted to see how it would look like:
    F(0) = std::cos(u[0] + u[1]) - 1 + 2 * u[0] - u[0];
    F(1) = std::sin(u[0] - u[1]) + 2 * u[1] - u[1];
    return 0;
  };

  FullMatrix<double> J_inverse(2, 2);

  kinsol.setup_jacobian = [&J_inverse](const VectorType &u,
                                       const VectorType &F) -> int {
    // We don't do any kind of set-up in this program, but we can at least
    // say that we're here
    deallog << "Setting up Jacobian system at u=(" << u[0] << ',' << u[1] << ")"
            << std::endl;

    FullMatrix<double> J(2, 2);
    J(0, 0) = -std::sin(u[0] + u[1]) + 2;
    J(0, 1) = -std::sin(u[0] + u[1]);
    J(1, 0) = std::cos(u[0] - u[1]);
    J(1, 1) = -std::cos(u[0] - u[1]) + 2;

    J_inverse.invert(J);

    return 0;
  };


  kinsol.solve_with_jacobian = [&J_inverse](const VectorType &rhs,
                                            VectorType &      dst,
                                            const double /*tolerance*/) -> int {
    deallog << "Solving Jacobian system with rhs=(" << rhs[0] << ',' << rhs[1]
            << ")" << std::endl;

    J_inverse.vmult(dst, rhs);

    return 0;
  };

  VectorType v(N);
  v(0) = 0.5;
  v(1) = 1.234;

  auto niter = kinsol.solve(v);
  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
