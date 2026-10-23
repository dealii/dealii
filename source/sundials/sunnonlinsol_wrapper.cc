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


#include <deal.II/base/config.h>

#include <deal.II/sundials/sunnonlinsol_wrapper.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SUNDIALS

namespace SUNDIALS
{
  namespace internal
  {
    SundialsNonlinearSolver
    make_fixed_point_nonlinear_solver(N_Vector y_template,
                                      int      anderson_acceleration_subspace,
                                      const InvocationContext &inv_ctx)
    {
#  if DEAL_II_SUNDIALS_VERSION_LT(6, 0, 0)
      (void)inv_ctx;
      return SundialsNonlinearSolver(
        SUNNonlinSol_FixedPoint(y_template, anderson_acceleration_subspace));
#  else
      return SundialsNonlinearSolver(
        SUNNonlinSol_FixedPoint(y_template,
                                anderson_acceleration_subspace,
                                inv_ctx.arkode_ctx));
#  endif
    }
  } // namespace internal
} // namespace SUNDIALS

#endif

DEAL_II_NAMESPACE_CLOSE
