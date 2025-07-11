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
 * Like the _04 test, but throw an error in the jacobian callbacks.
 * Unlike for the SNES case in petsc_snes_07.cc, recoverable errors
 * are fully recoverable and lead to a smaller step size selected.
 * Here we also test the case of mixing the explicit and implicit interface.
 * That's why we first run the case without raising any exception.
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
                   bool                                           _testerr,
                   bool                                           _testerrrec,
                   int                                            _testcase)
    : time_stepper(data)
    , kappa(_kappa)
    , last_eval_time(0)
    , testerr(_testerr)
    , testerrrecoverable(_testerrrec)
    , testcase(_testcase)
  {
    throw_error = false;

    time_stepper.explicit_function =
      [&](const real_type t, const VectorType &y, VectorType &res) -> void {
      deallog << "Evaluating explicit function at t=" << t << std::endl;

      if (testerr && (t > last_eval_time + 0.1))
        {
          deallog << "Time step too large: last_eval_time=" << last_eval_time
                  << ", t=" << t << std::endl;
          testerr     = false;
          throw_error = true;
        }

      res(0) = -kappa * y(0);
      res.compress(VectorOperation::insert);

      last_eval_time = t;
    };

    time_stepper.monitor = [&](const real_type   t,
                               const VectorType &sol,
                               const unsigned int) -> void {
      deallog << "Intermediate output:" << std::endl;
      deallog << "  t =" << t << std::endl;
      deallog << "  y =" << sol[0] << "  (exact: " << std::exp(-kappa * t)
              << ')' << std::endl;
    };

    if (!testcase)
      {
        time_stepper.implicit_jacobian = [&](const real_type   t,
                                             const VectorType &y,
                                             const VectorType &y_dot,
                                             const real_type   shift,
                                             MatrixType       &A,
                                             MatrixType       &P) -> void {
          deallog << "Evaluating implicit Jacobian at t=" << t << std::endl;
          if (throw_error)
            {
              throw_error = false;
              if (testerrrecoverable)
                {
                  deallog
                    << "Throwing a recoverable exception from implicit_jacobian."
                    << std::endl;
                  throw RecoverableUserCallbackError();
                }
              else
                {
                  deallog
                    << "Throwing a unrecoverable exception from implicit_jacobian."
                    << std::endl;
                  throw std::runtime_error(
                    "Unrecoverable error in implicit_jacobian.");
                }
            }
          P.set(0, 0, shift + kappa);
          P.compress(VectorOperation::insert);
        };
      }
    else
      {
        time_stepper.setup_jacobian = [&](const real_type   t,
                                          const VectorType &y,
                                          const VectorType &y_dot,
                                          const real_type   shift) -> void {
          deallog << "Setting up Jacobian at t=" << t << std::endl;
          myshift = shift;
          if (throw_error && testcase == 1)
            {
              throw_error = false;
              if (testerrrecoverable)
                {
                  deallog
                    << "Throwing a recoverable exception from setup_jacobian."
                    << std::endl;
                  throw RecoverableUserCallbackError();
                }
              else
                {
                  deallog
                    << "Throwing a unrecoverable exception from setup_jacobian."
                    << std::endl;
                  throw std::runtime_error(
                    "Unrecoverable error in setup_jacobian.");
                }
            }
        };

        time_stepper.solve_with_jacobian = [&](const VectorType &src,
                                               VectorType       &dst) -> void {
          deallog << "Solving with Jacobian" << std::endl;
          if (throw_error && testcase == 2)
            {
              throw_error = false;
              if (testerrrecoverable)
                {
                  deallog
                    << "Throwing a recoverable exception from solve_with_jacobian."
                    << std::endl;
                  throw RecoverableUserCallbackError();
                }
              else
                {
                  deallog
                    << "Throwing a unrecoverable exception from solve_with_jacobian."
                    << std::endl;
                  throw std::runtime_error(
                    "Unrecoverable error in setup_jacobian.");
                }
            }
          dst(0) = src(0) / (myshift + kappa);
          dst.compress(VectorOperation::insert);
        };
      }
  }

  void
  run()
  {
    auto t = time_stepper.get_time();

    VectorType y(MPI_COMM_SELF, 1, 1);
    y[0] = 1;
    y.compress(VectorOperation::insert);

    try
      {
        auto nt = time_stepper.solve(y);
        deallog << "# Number of steps taken: " << nt << std::endl;
      }
    catch (const std::exception &exc)
      {
        deallog << "Time stepper aborted with an exception:" << std::endl
                << exc.what() << std::endl;
      }
  }

private:
  TimeStepper time_stepper;
  real_type   kappa;
  real_type   myshift;
  double      last_eval_time;
  bool        testerr;
  bool        testerrrecoverable;
  int         testcase;
  bool        throw_error;
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

  // This test triggers false positives in FPE trapping for some versions of
  // PETSc
#if DEAL_II_PETSC_VERSION_LT(3, 19, 2) && defined(DEBUG) && \
  defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  PetscErrorCode ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);
  (void)ierr;
#endif

  for (int seterri = 0; seterri < 2; seterri++)
    {
      bool seterr = seterri ? true : false;

      for (int testcase = 0; testcase < 3; testcase++)
        {
          deallog << "# Testcase " << testcase;
          if (seterr)
            deallog << " with recoverable error";
          deallog << std::endl;
          ExponentialDecay ode_expl(1.0, data, seterr, true, testcase);
          ode_expl.run();
          if (seterr)
            {
              deallog << "# Testcase " << testcase
                      << " with unrecoverable error " << std::endl;
              ExponentialDecay ode_expl(1.0, data, seterr, false, testcase);
              ode_expl.run();
            }
        }
    }
}
