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
 * Solves the equations of exponential decay:
 *
 * u' = -k u
 * u (t0) = 1
 *
 * using the PETScWrappers::TimeStepper class
 * that interfaces PETSc TS ODE solver object.
 * The ODE can be an EXPLICIT first order ode:
 *
 * y[0]' = - k   y[0]     -> y' = G(y,t)
 *
 * or specified in IMPLICIT form:
 *
 * y[0]' +  k  y[1]  = 0  -> F(y',y,t) = 0
 *
 * The exact solution is, in either case,
 *
 * y[0](t) = exp(-kt)
 *
 * We use the same class to test both formulations.
 * The class also supports a third approach where
 * control for linear systems generation and
 * solution is completely on users. In this case
 * users are in charge of solving for the
 * implicit Jacobian.
 * Here we also test the functionalities of the
 * various user callbacks available.
 * The model used to wrap PETSc's TS is the same
 * used for the nonlinear solver SNES. Check
 * petsc_snes_00.cc for additional information.
 *
 */
using VectorType  = PETScWrappers::MPI::Vector;
using MatrixType  = PETScWrappers::MatrixBase;
using TimeStepper = PETScWrappers::TimeStepper<VectorType, MatrixType>;
using real_type   = TimeStepper::real_type;

class ExponentialDecay
{
public:
  ExponentialDecay(real_type                                      _kappa,
                   const typename PETScWrappers::TimeStepperData &data,
                   bool                                           setjac,
                   bool                                           implicit,
                   bool                                           user)
    : time_stepper(data)
    , kappa(_kappa)
  {
    // In this case we use the implicit form
    if (implicit)
      {
        time_stepper.implicit_function = [&](const real_type   t,
                                             const VectorType &y,
                                             const VectorType &y_dot,
                                             VectorType       &res) -> void {
          deallog << "Evaluating implicit function at t=" << t << std::endl;
          res(0) = y_dot(0) + kappa * y(0);
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
              deallog << "Evaluating implicit Jacobian at t=" << t << std::endl;
              P.set(0, 0, shift + kappa);
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
              deallog << "Setting up Jacobian at t=" << t << std::endl;
              myshift = shift;
            };

            // In the solve phase we se the stored shift to solve
            // for the implicit Jacobian system
            time_stepper.solve_with_jacobian = [&](const VectorType &src,
                                                   VectorType &dst) -> void {
              deallog << "Solving with Jacobian" << std::endl;
              dst(0) = src(0) / (myshift + kappa);
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
          deallog << "Evaluating explicit function at t=" << t << std::endl;
          res(0) = -kappa * y(0);
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
              deallog << "Evaluating explicit Jacobian at t=" << t << std::endl;
              P.set(0, 0, -kappa);
              P.compress(VectorOperation::insert);
            };
          }
      }

    // Monitoring routine. Here we print diagnostic for the exact
    // solution to the log file.
    time_stepper.monitor = [&](const real_type   t,
                               const VectorType &sol,
                               const unsigned int /*step_number*/) -> void {
      deallog << "Intermediate output:" << std::endl;
      deallog << "  t =" << t << std::endl;
      deallog << "  y =" << sol[0] << "  (exact: " << std::exp(-kappa * t)
              << ')' << std::endl;
    };

    // This callback is invoked after a successful stage.
    // Here we only print that the callback is invoked.
    time_stepper.update_constrained_components = [&](const real_type t,
                                                     VectorType &) -> void {
      deallog << "Distribute at time " << t << std::endl;
    };

    // This callback is used to decide to remesh.
    time_stepper.decide_and_prepare_for_remeshing =
      [&](const real_type    t,
          const unsigned int step,
          const VectorType &) -> bool {
      deallog << "Prepare at time " << t << " and step " << step << std::endl;
      return (step && step % 5 == 0);
    };

    // This callback is called if decide_and_prepare_for_remeshing returns true.
    time_stepper.transfer_solution_vectors_to_new_mesh =
      [&](const double /* t */,
          const std::vector<VectorType> &all_in,
          std::vector<VectorType>       &all_out) -> void {
      deallog << "Interpolate" << std::endl;
      for (auto &v : all_in)
        all_out.push_back(v);
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

    VectorType y(MPI_COMM_SELF, 1, 1);
    y[0] = 1;
    y.compress(VectorOperation::insert);

    // Integrate the ODE.
    auto nt = time_stepper.solve(y);
    deallog << "# Number of steps taken: " << nt << std::endl;
  }

private:
  TimeStepper time_stepper;
  real_type   kappa;   // Defines the decay
  real_type   myshift; // Used by the user solve
};


int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  PETScWrappers::TimeStepperData data;
  ParameterHandler               prm;

  data.add_parameters(prm);
  deallog << "# Default Parameters" << std::endl;
  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::ShortPRM);

  std::ifstream ifile(SOURCE_DIR "/petsc_ts_03_in.prm");
  prm.parse_input(ifile);

  deallog << "# Testing Parameters" << std::endl;
  prm.print_parameters(deallog.get_file_stream(), ParameterHandler::ShortPRM);

  for (int setjaci = 0; setjaci < 2; setjaci++)
    {
      bool setjac = setjaci ? true : false;

      {
        deallog << "# Test explicit interface (J " << setjac << ")"
                << std::endl;
        ExponentialDecay ode_expl(1.0, data, setjac, false, false);
        ode_expl.run();
      }

      {
        deallog << "# Test implicit interface (J " << setjac << ")"
                << std::endl;
        ExponentialDecay ode_impl(1.0, data, setjac, true, false);
        ode_impl.run();
      }
    }

  {
    deallog << "# Test user interface" << std::endl;
    ExponentialDecay ode_user(1.0, data, true, true, true);
    ode_user.run();
  }
}
