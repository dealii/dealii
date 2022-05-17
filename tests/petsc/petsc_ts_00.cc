//-----------------------------------------------------------
//
//    Copyright (C) 2022 by the deal.II authors
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
//
// Author: Stefano Zampini, King Abdullah University of Science and Technology.

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/petsc_matrix_base.h>
#include <deal.II/lac/petsc_ts.h>
#include <deal.II/lac/petsc_vector.h>

#include <cmath>

#include "../tests.h"

/**
 * Solves the harmonic oscillator problem
 *
 * u'' = -k^2 u
 * u (t0) = sin(k * t0)
 * u'(t0) = k cos (k* t0)
 *
 * using the PETScWrappers::TimeStepper class
 * that interfaces PETSc TS ODE solver object.
 * The ODE can be an EXPLICIT first order ode:
 *
 * y[0]' =       y[1]   \
 *                       -> y' = G(y,t)
 * y[1]' = - k^2 y[0]   /
 *
 * or specified in IMPLICIT form:
 *
 * y[0]' -     y[1]  = 0 \
 *                        -> F(y',y,t) = 0
 * y[1]' + k^2 y[0]  = 0 /
 *
 * The exact solution is
 *
 * y[0](t) = sin(k t)
 * y[1](t) = k cos(k t)
 *
 * We use the same class to test both formulations.
 * The class also supports a third approach where
 * control for linear systems generation and
 * solution is completely on users. In this case
 * users are in charge of solving for the
 * implicit Jacobian.
 * The model used to wrap PETSc's TS is the same
 * used for the nonlinear solver SNES. Check
 * petsc_snes_00.cc for additional information.
 *
 */
using VectorType  = PETScWrappers::MPI::Vector;
using MatrixType  = PETScWrappers::MatrixBase;
using TimeStepper = PETScWrappers::TimeStepper<VectorType, MatrixType>;

class HarmonicOscillator
{
public:
  HarmonicOscillator(double                                         _kappa,
                     const typename PETScWrappers::TimeStepperData &data,
                     bool                                           implicit,
                     bool                                           user,
                     std::ostream &                                 _out)
    : time_stepper(data)
    , out(_out)
    , kappa(_kappa)
  {
    // In this case we use the implicit form
    if (implicit)
      {
        time_stepper.implicit_function = [&](const double      t,
                                             const VectorType &y,
                                             const VectorType &y_dot,
                                             VectorType &      res) -> int {
          res(0) = y_dot(0) - y(1);
          res(1) = y_dot(1) + kappa * kappa * y(0);
          res.compress(VectorOperation::insert);
          return 0;
        };

        // We either have the possibility of using PETSc standard
        // functionalities for Jacobian setup and solve, or take full control of
        // the Jacobian solutions. The latter aligns more with deal.II tutorials
        // For our testing class, 'user' discriminate between those two
        // approaches.
        if (!user)
          {
            // Callback for Jacobian evaluation
            // This callback populates the P matrix with
            //     shift * dF/dydot + dF/dy
            // A and P are the same matrix if not otherwise specified
            // PETSc can compute a Jacobian by finite-differencing the residual
            // evaluation, and it can be selected at command line.
            time_stepper.implicit_jacobian = [&](const double      t,
                                                 const VectorType &y,
                                                 const VectorType &y_dot,
                                                 const double      shift,
                                                 MatrixType &      A,
                                                 MatrixType &      P) -> int {
              P.set(0, 0, shift);
              P.set(0, 1, -1);
              P.set(1, 0, kappa * kappa);
              P.set(1, 1, shift);
              P.compress(VectorOperation::insert);
              return 0;
            };
          }
        else
          {
            // In this code block, users are completely in charge of
            // setting up the Jacobian system and solve for it.
            // In this example we only store the solver shift
            // during setup.
            time_stepper.setup_jacobian = [&](const double      t,
                                              const VectorType &y,
                                              const VectorType &y_dot,
                                              const double      shift) -> int {
              myshift = shift;
              return 0;
            };

            // In the solve phase we se the stored shift to solve
            // for the implicit Jacobian system
            time_stepper.solve_for_jacobian_system =
              [&](const VectorType &src, VectorType &dst) -> int {
              auto sf = 1. / (kappa * kappa + myshift * myshift);
              dst(0)  = sf * (myshift * src(0) + src(1));
              dst(1)  = sf * (-kappa * kappa * src(0) + myshift * src(1));
              dst.compress(VectorOperation::insert);
              return 0;
            };
          }
      }
    else
      { // Here we instead use the explicit form
        // This is the only function one would populate in case an explicit
        // solver is used.
        time_stepper.explicit_function =
          [&](const double t, const VectorType &y, VectorType &res) -> int {
          res(0) = y(1);
          res(1) = -kappa * kappa * y(0);
          res.compress(VectorOperation::insert);
          return 0;
        };

        // The explicit Jacobian callback is not needed in case
        // an explicit solver is used. We need it when moving to
        // an implicit solver like in the parameter file used for
        // this test
        time_stepper.explicit_jacobian = [&](const double      t,
                                             const VectorType &y,
                                             MatrixType &      A,
                                             MatrixType &      P) -> int {
          P.set(0, 0, 0);
          P.set(0, 1, 1);
          P.set(1, 0, -kappa * kappa);
          P.set(1, 1, 0);
          P.compress(VectorOperation::insert);
          return 0;
        };
      }

    // Monitoring routine. Here we print diagnostic for the exact
    // solution to the log file.
    time_stepper.monitor = [&](const double       t,
                               const VectorType & y,
                               const unsigned int step_number) -> int {
      std::vector<double> exact(2);
      exact[0] = std::sin(kappa * t);
      exact[1] = kappa * std::cos(kappa * t);
      out << t << " " << y(0) << " (" << exact[0] << ")"
          << " " << y(1) << " (" << exact[1] << ")" << std::endl;
      return 0;
    };
  }

  void
  run()
  {
    // Set initial conditions
    // This test uses a nonzero initial time
    // We can access the time value with the
    // get_time method
    auto t = time_stepper.get_time();

    VectorType y(MPI_COMM_SELF, 2, 2);
    y[0] = std::sin(kappa * t);
    y[1] = kappa * std::cos(kappa * t);
    y.compress(VectorOperation::insert);

    // Integrate the ODE.
    // Here we do not pass any matrix and our runs will use the default
    // matrix type in PETSc (MATDENSE) in case a PETSc interface is used.
    time_stepper.solve(y);
  }

private:
  TimeStepper   time_stepper;
  std::ostream &out;     // Used by the monitoring routine
  double        kappa;   // Defines the oscillator
  double        myshift; // Used by the user solve
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream out("output");

  PETScWrappers::TimeStepperData data;
  ParameterHandler               prm;

  data.add_parameters(prm);
  out << "# Default Parameters" << std::endl;
  prm.print_parameters(out, ParameterHandler::ShortText);

  std::ifstream ifile(SOURCE_DIR "/petsc_ts_00_in.prm");
  prm.parse_input(ifile);

  out << "# Testing Parameters" << std::endl;
  prm.print_parameters(out, ParameterHandler::ShortText);

  {
    out << "# Test explicit interface" << std::endl;
    HarmonicOscillator ode_expl(1.0, data, false, false, out);
    ode_expl.run();
  }

  {
    out << "# Test implicit interface" << std::endl;
    HarmonicOscillator ode_impl(1.0, data, true, false, out);
    ode_impl.run();
  }

  {
    out << "# Test user interface" << std::endl;
    HarmonicOscillator ode_user(1.0, data, true, true, out);
    ode_user.run();
  }
  return 0;
}
