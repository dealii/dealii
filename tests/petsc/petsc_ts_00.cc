// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
using real_type   = TimeStepper::real_type;

class HarmonicOscillator
{
public:
  HarmonicOscillator(real_type                                      _kappa,
                     const typename PETScWrappers::TimeStepperData &data,
                     bool                                           setjac,
                     bool                                           implicit,
                     bool                                           user,
                     std::ostream                                  &_out)
    : time_stepper(data)
    , out(_out)
    , kappa(_kappa)
  {
    // In this case we use the implicit form
    if (implicit)
      {
        time_stepper.implicit_function = [&](const real_type   t,
                                             const VectorType &y,
                                             const VectorType &y_dot,
                                             VectorType       &res) -> void {
          res(0) = y_dot(0) - y(1);
          res(1) = y_dot(1) + kappa * kappa * y(0);
          res.compress(VectorOperation::insert);
        };

        // We either have the possibility of using PETSc standard
        // functionalities for Jacobian setup and solve, or take full control of
        // the Jacobian solutions. The latter aligns more with deal.II tutorials
        // For our testing class, 'user' discriminate between those two
        // approaches.
        if (!user && setjac)
          {
            // Callback for Jacobian evaluation
            // This callback populates the P matrix with
            //     shift * dF/dydot + dF/dy
            // A and P are the same matrix if not otherwise specified
            // PETSc can compute a Jacobian by finite-differencing the residual
            // evaluation, and it can be selected at command line.
            time_stepper.implicit_jacobian = [&](const real_type   t,
                                                 const VectorType &y,
                                                 const VectorType &y_dot,
                                                 const real_type   shift,
                                                 MatrixType       &A,
                                                 MatrixType       &P) -> void {
              P.set(0, 0, shift);
              P.set(0, 1, -1);
              P.set(1, 0, kappa * kappa);
              P.set(1, 1, shift);
              P.compress(VectorOperation::insert);
            };
          }
        else if (user)
          {
            // In this code block, users are completely in charge of
            // setting up the Jacobian system and solve for it.
            // In this example we only store the solver shift
            // during setup.
            time_stepper.setup_jacobian = [&](const real_type   t,
                                              const VectorType &y,
                                              const VectorType &y_dot,
                                              const real_type   shift) -> void {
              myshift = shift;
            };

            // In the solve phase we se the stored shift to solve
            // for the implicit Jacobian system
            time_stepper.solve_with_jacobian = [&](const VectorType &src,
                                                   VectorType &dst) -> void {
              auto sf = 1. / (kappa * kappa + myshift * myshift);
              dst(0)  = sf * (myshift * src(0) + src(1));
              dst(1)  = sf * (-kappa * kappa * src(0) + myshift * src(1));
              dst.compress(VectorOperation::insert);
            };
          }
      }
    else
      { // Here we instead use the explicit form
        // This is the only function one would populate in case an explicit
        // solver is used.
        time_stepper.explicit_function =
          [&](const real_type t, const VectorType &y, VectorType &res) -> void {
          res(0) = y(1);
          res(1) = -kappa * kappa * y(0);
          res.compress(VectorOperation::insert);
        };

        // The explicit Jacobian callback is not needed in case
        // an explicit solver is used. We need it when moving to
        // an implicit solver like in the parameter file used for
        // this test
        if (setjac)
          {
            time_stepper.explicit_jacobian = [&](const real_type   t,
                                                 const VectorType &y,
                                                 MatrixType       &A,
                                                 MatrixType       &P) -> void {
              P.set(0, 0, 0);
              P.set(0, 1, 1);
              P.set(1, 0, -kappa * kappa);
              P.set(1, 1, 0);
              P.compress(VectorOperation::insert);
            };
          }
      }

    // Monitoring routine. Here we print diagnostic for the exact
    // solution to the log file.
    time_stepper.monitor = [&](const real_type    t,
                               const VectorType  &y,
                               const unsigned int step_number) -> void {
      std::vector<real_type> exact(2);
      exact[0] = std::sin(kappa * t);
      exact[1] = kappa * std::cos(kappa * t);
      out << t << " " << y(0) << " (" << exact[0] << ")"
          << " " << y(1) << " (" << exact[1] << ")" << std::endl;
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
    auto nt = time_stepper.solve(y);
    out << "# Number of steps taken: " << nt << std::endl;
  }

private:
  TimeStepper   time_stepper;
  std::ostream &out;     // Used by the monitoring routine
  real_type     kappa;   // Defines the oscillator
  real_type     myshift; // Used by the user solve
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
  prm.print_parameters(out, ParameterHandler::ShortPRM);

  std::ifstream ifile(SOURCE_DIR "/petsc_ts_00_in.prm");
  prm.parse_input(ifile);

  out << "# Testing Parameters" << std::endl;
  prm.print_parameters(out, ParameterHandler::ShortPRM);

  for (int setjaci = 0; setjaci < 2; setjaci++)
    {
      bool setjac = setjaci ? true : false;

      {
        out << "# Test explicit interface (J " << setjac << ")" << std::endl;
        HarmonicOscillator ode_expl(1.0, data, setjac, false, false, out);
        ode_expl.run();
      }

      {
        out << "# Test implicit interface (J " << setjac << ")" << std::endl;
        HarmonicOscillator ode_impl(1.0, data, setjac, true, false, out);
        ode_impl.run();
      }
    }

  {
    out << "# Test user interface" << std::endl;
    HarmonicOscillator ode_user(1.0, data, true, true, true, out);
    ode_user.run();
  }
}
