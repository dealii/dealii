//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2023 by the deal.II authors
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

// Test explicit time stepper where we throw a recoverable error from
// the function that computes the right hand side. The test case
// solves the problem y'(t)=-y(t) with exact solution
// y(t)=exp(-t). The right hand side function throws an exception if
// evaluated with y<0, and we make sure that that happens by starting
// with an overly large initial time step.
int
main()
{
  initlog();

  using VectorType = Vector<double>;

  ParameterHandler                             prm;
  SUNDIALS::ARKode<VectorType>::AdditionalData data;
  data.add_parameters(prm);

  // Set to true to reset input file.
  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_09_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortText);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_09_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKode<VectorType> ode(data);

  double kappa = 1.0;

  ode.explicit_function =
    [&](const double t, const VectorType &y, VectorType &ydot) {
      deallog << "Evaluating right hand side at t=" << t << " with y=" << y[0]
              << std::endl;

      // Error out in a recoverable way if asked to evaluate at a
      // point where y<0. This can happen with explicit methods if the
      // time step is too large.
      if (y[0] < 0)
        {
          deallog << "Reporting a recoverable error." << std::endl;
          throw RecoverableUserCallbackError();
        }

      ydot[0] = -y[0];
    };

  ode.output_step =
    [&](const double t, const VectorType &sol, const unsigned int step_number) {
      deallog << "Monitor called at t=" << t << ' ' << sol[0] << std::endl;
    };

  Vector<double> y(1);
  y[0] = 1;

  try
    {
      ode.solve_ode(y);
    }
  catch (const std::exception &exc)
    {
      deallog << "ARKODE did not succeed with the following error message:"
              << std::endl
              << exc.what() << std::endl;
    }
}
