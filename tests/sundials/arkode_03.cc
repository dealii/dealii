//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2018 by the deal.II authors
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
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/sundials/arkode.h>

#include "../tests.h"


// Test implicit-explicit time stepper, no jacobian. Use L:d:V (in serial)
// Brusselator benchmark

/**
 * This test problem is called "brusselator", and is a typical benchmark for
 * ODE solvers. This problem has 3 dependent variables u, v and w, that depend
 * on the independent variable t via the IVP system
 *
 * du/dt = a − (w + 1)u + v u^2
 * dv/dt = w u − v u^2
 * dw/dt = (b − w)/eps -w u
 *
 * We integrate over the interval 0 ≤ t ≤ 10, with the initial conditions
 *
 * u(0) = 3.9, v(0) = 1.1, w(0) = 2.8,
 *
 * and parameters
 *
 * a = 1.2, b = 2.5, and eps = 10−5
 *
 * The implicit part only contains the stiff part of the problem (the part with
 * eps in right hand side of the third equation).
 */
int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  using VectorType = LinearAlgebra::distributed::Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_03.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortText);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_03.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  // Parameters
  double u0 = 3.9, v0 = 1.1, w0 = 2.8, a = 1.2, b = 2.5, eps = 1e-5;

  ode.implicit_function =
    [&](double, const VectorType &y, VectorType &ydot) -> int {
    ydot[0] = 0;
    ydot[1] = 0;
    ydot[2] = (b - y[2]) / eps;
    return 0;
  };


  ode.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) -> int {
    ydot[0] = a - (y[2] + 1) * y[0] + y[1] * y[0] * y[0];
    ydot[1] = y[2] * y[0] - y[1] * y[0] * y[0];
    ydot[2] = -y[2] * y[0];
    return 0;
  };

  ode.output_step = [&](const double       t,
                        const VectorType & sol,
                        const unsigned int step_number) -> int {
    // limit the output to every 10th step and increase the precision to make
    // the test more robust
    if (step_number % 10 == 0)
      deallog << t << " " << std::setprecision(10) << sol[0] << " " << sol[1]
              << " " << sol[2] << std::endl;
    return 0;
  };

  VectorType y(3);
  y[0] = u0;
  y[1] = v0;
  y[2] = w0;
  ode.solve_ode(y);
  return 0;
}
