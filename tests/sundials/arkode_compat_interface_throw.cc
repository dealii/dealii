// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/arkode.h>
#include <deal.II/sundials/arkode_stepper.h>

#include "../tests.h"

// Verify that every public FunctionProxy member of SUNDIALS::ARKode throws
// when ARKode is constructed with an explicit ARKodeStepper argument.
//
// When ARKode is constructed with a user-supplied stepper the proxy target
// pointers are nullptr, because the callbacks live on the stepper, not on
// ARKode.  Assigning a std::function to a proxy whose target is nullptr must
// trigger an AssertThrow (ExcMessage) in both debug and release builds.
// Each sub-test below checks one proxy member, and the test is exhaustive
// over all twelve public FunctionProxy fields.

using VectorType = Vector<double>;

// Helper: try assigning to a proxy and report the result.
#define CHECK_PROXY(name, assignment)                                  \
  do                                                                   \
    {                                                                  \
      try                                                              \
        {                                                              \
          assignment;                                                  \
          deallog << #name ": ERROR – expected exception not thrown" \
                  << std::endl;                                        \
        }                                                              \
      catch (const ExceptionBase &)                                    \
        {                                                              \
          deallog << #name ": OK" << std::endl;                        \
        }                                                              \
    }                                                                  \
  while (false)

int
main()
{
  initlog();

  SUNDIALS::ARKStepper<VectorType> stepper;
  SUNDIALS::ARKode<VectorType>     ode(stepper);

  CHECK_PROXY(
    explicit_function,
    ode.explicit_function = [](const double, const VectorType &, VectorType &) {
    });

  CHECK_PROXY(
    implicit_function,
    ode.implicit_function = [](const double, const VectorType &, VectorType &) {
    });

  CHECK_PROXY(
    mass_times_vector,
    ode.mass_times_vector = [](const double, const VectorType &, VectorType &) {
    });

  CHECK_PROXY(
    mass_times_setup, ode.mass_times_setup = [](const double) {});

  CHECK_PROXY(
    jacobian_times_vector,
    ode.jacobian_times_vector = [](const VectorType &,
                                   VectorType &,
                                   const double,
                                   const VectorType &,
                                   const VectorType &) {});

  CHECK_PROXY(
    jacobian_times_setup,
    ode.jacobian_times_setup =
      [](const double, const VectorType &, const VectorType &) {});

  CHECK_PROXY(
    solve_linearized_system,
    ode.solve_linearized_system =
      [](SUNDIALS::SundialsOperator<VectorType> &,
         SUNDIALS::SundialsPreconditioner<VectorType> &,
         VectorType &,
         const VectorType &,
         double) {});

  CHECK_PROXY(
    solve_mass,
    ode.solve_mass = [](SUNDIALS::SundialsOperator<VectorType> &,
                        SUNDIALS::SundialsPreconditioner<VectorType> &,
                        VectorType &,
                        const VectorType &,
                        double) {});

  CHECK_PROXY(
    jacobian_preconditioner_solve,
    ode.jacobian_preconditioner_solve = [](const double,
                                           const VectorType &,
                                           const VectorType &,
                                           const VectorType &,
                                           VectorType &,
                                           const double,
                                           const double,
                                           const int) {});

  CHECK_PROXY(
    jacobian_preconditioner_setup,
    ode.jacobian_preconditioner_setup = [](const double,
                                           const VectorType &,
                                           const VectorType &,
                                           const int,
                                           int &,
                                           const double) {});

  CHECK_PROXY(
    mass_preconditioner_solve,
    ode.mass_preconditioner_solve = [](const double,
                                       const VectorType &,
                                       VectorType &,
                                       const double,
                                       const int) {});

  CHECK_PROXY(
    mass_preconditioner_setup,
    ode.mass_preconditioner_setup = [](const double) {});
}
