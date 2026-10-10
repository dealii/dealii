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

#include "arkode_restart_common.h"

// Test ARKode::solver_should_restart() with ARKStepper, with an explicit and
// with an implicit right hand side. The implicit one is solved with a user
// provided linear solver, with SUNDIALS' own linear solver and with SUNDIALS'
// fixed-point solver. A restart has to release all of them.

int
main()
{
  initlog();

  using Stepper = SUNDIALS::ARKStepper<VectorType>;

  const auto set_jacobian = [](Stepper &stepper) {
    stepper.implicit_function     = rhs;
    stepper.jacobian_times_vector = [](const VectorType &v,
                                       VectorType       &Jv,
                                       double,
                                       const VectorType &,
                                       const VectorType &) {
      Jv = v;
      Jv *= -1.;
    };
  };

  run_restart_test<Stepper>("ARKStepper explicit");

  run_restart_test<Stepper>("ARKStepper implicit", [&](Stepper &stepper) {
    set_jacobian(stepper);
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

  run_restart_test<Stepper>("ARKStepper implicit with SUNDIALS linear solver",
                            set_jacobian);

  run_restart_test<Stepper>("ARKStepper implicit with fixed-point solver",
                            [](Stepper &stepper) {
                              stepper.implicit_function = rhs;
                            });
}
