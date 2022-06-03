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

// Solve a nonlinear system using fixed point iteration, and Anderson
// acceleration

/**
 * The following is a simple example problem, with the coding
 * needed for its solution by the accelerated fixed point solver in
 * KINSOL.
 * The problem is from chemical kinetics, and consists of solving
 * the first time step in a Backward Euler solution for the
 * following three rate equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2
 *    dy3/dt = 3.e2*(y2)^2
 * on the interval from t = 0.0 to t = 0.1, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * Run statistics (optional outputs) are printed at the end.
 */
int
main()
{
  initlog();

  using VectorType = Vector<double>;

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/kinsol_fixed_point.prm");
  prm.parse_input(ifile);

  // Size of the problem
  unsigned int N = 3;

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  // Robert example
  kinsol.iteration_function = [](const VectorType &u, VectorType &F) -> int {
    const double dstep = 0.1;
    const double y10   = 1.0;
    const double y20   = 0.0;
    const double y30   = 0.0;

    const double yd1 = dstep * (-0.04 * u[0] + 1.0e4 * u[1] * u[2]);
    const double yd3 = dstep * 3.0e2 * u[1] * u[1];

    F[0] = yd1 + y10;
    F[1] = -yd1 - yd3 + y20;
    F[2] = yd3 + y30;

    return 0;
  };

  VectorType v(N);
  v[0]       = 1;
  auto niter = kinsol.solve(v);
  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
