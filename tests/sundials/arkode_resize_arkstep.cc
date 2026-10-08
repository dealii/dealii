// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include "arkode_resize_common.h"

// Test ARKode::solver_should_resize() with ARKStepper, with an explicit and
// with an implicit right hand side. The latter uses a linear solver that has
// to be re-attached to the resized vectors.

int
main()
{
  initlog();

  using Stepper = SUNDIALS::ARKStepper<VectorType>;

  run_resize_test<Stepper>("ARKStepper explicit");

  run_resize_test<Stepper>("ARKStepper implicit", [](Stepper &stepper) {
    stepper.implicit_function     = rhs;
    stepper.jacobian_times_vector = [](const VectorType &v,
                                       VectorType       &Jv,
                                       double,
                                       const VectorType &,
                                       const VectorType &) {
      Jv = v;
      Jv *= -1.;
    };
    stepper.solve_linearized_system =
      [](SUNDIALS::SundialsOperator<VectorType> &op,
         SUNDIALS::SundialsPreconditioner<VectorType> &,
         VectorType       &x,
         const VectorType &b,
         double            tol) {
        SolverControl        control(100, tol);
        SolverCG<VectorType> solver_cg(control);
        solver_cg.solve(op, x, b, PreconditionIdentity());
      };
  });
}
