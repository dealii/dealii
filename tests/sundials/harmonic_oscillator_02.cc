//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
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

#include <deal.II/sundials/arkode.h>

#include "../tests.h"

// Test explicit time stepper. Only implements explicit_function.

/**
 * Solve the Harmonic oscillator problem.
 *
 * u'' = -k^2 u
 * u (0) = 0
 * u'(0) = k
 *
 * write in terms of a first order ode:
 *
 * y[0]' =       y[1]
 * y[1]' = - k^2 y[0]
 *
 * That is
 *
 * y' = A y
 *
 * A = [ 0 , 1; -k^2, 0 ]
 *
 * y_0  = 0, k
 *
 * The exact solution is
 *
 * y[0](t) = sin(k t)
 * y[1](t) = k cos(k t)
 *
 */
int
main(int argc, char **argv)
{
  std::ofstream out("output");

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  typedef Vector<double> VectorType;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  // Set to true to reset input file.
  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/harmonic_oscillator_02.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortText);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/harmonic_oscillator_02.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  ode.reinit_vector = [&](VectorType &v) { v.reinit(2); };

  double kappa = 1.0;

  ode.explicit_function =
    [&](double, const VectorType &y, VectorType &ydot) -> int {
    ydot[0] = y[1];
    ydot[1] = -kappa * kappa * y[0];
    return 0;
  };

  ode.output_step = [&](const double       t,
                        const VectorType & sol,
                        const unsigned int step_number) -> int {
    out << t << " " << sol[0] << " " << sol[1] << std::endl;
    return 0;
  };

  Vector<double> y(2);
  y[0] = 0;
  y[1] = kappa;
  ode.solve_ode(y);
  return 0;
}
