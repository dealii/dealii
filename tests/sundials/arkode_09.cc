// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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

  SUNDIALS::ARKStepper<VectorType>::AdditionalData stepper_data;
  stepper_data.add_parameters(prm);

  // Set to true to reset input file.
  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/arkode_09_in.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortPRM);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/arkode_09_in.prm");
  prm.parse_input(ifile);

  SUNDIALS::ARKStepper<VectorType> stepper(stepper_data);
  SUNDIALS::ARKode<VectorType>     ode(stepper, data);

  double kappa = 1.0;

  bool explicit_function_called = false;

  stepper.explicit_function =
    [&](const double t, const VectorType &y, VectorType &ydot) {
      if (!explicit_function_called)
        {
          deallog << "Right hand side callback has been called." << std::endl;
          explicit_function_called = true;
        }

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
