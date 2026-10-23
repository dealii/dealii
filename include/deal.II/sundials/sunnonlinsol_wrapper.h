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

#ifndef dealii_sundials_sunnonlinsol_wrapper_h
#define dealii_sundials_sunnonlinsol_wrapper_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/sundials/invocation_context.h>
#  include <deal.II/sundials/sundials_handle.h>

#  include <sundials/sundials_nonlinearsolver.h>
#  include <sundials/sundials_nvector.h>

#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_SUNDIALS
namespace SUNDIALS
{
  namespace internal
  {
    /**
     * Traits for SundialsHandle that free a SUNNonlinearSolver with
     * SUNNonlinSolFree().
     */
    struct SundialsNonlinearSolverTraits
    {
      using Handle = SUNNonlinearSolver;

      static void
      free(Handle solver)
      {
        SUNNonlinSolFree(solver);
      }
    };

    /**
     * Owning handle for any SUNNonlinearSolver created by SUNDIALS, for
     * example by SUNNonlinSol_FixedPoint() or SUNNonlinSol_Newton().
     */
    using SundialsNonlinearSolver =
      SundialsHandle<SundialsNonlinearSolverTraits>;

    /**
     * Create SUNDIALS' fixed-point nonlinear solver.
     *
     * @param y_template A vector that is used as a template to create the
     *   internal vectors of the solver. It is not needed after this function
     *   returns.
     * @param anderson_acceleration_subspace Size of the Anderson acceleration
     *   subspace.
     * @param inv_ctx The invocation context that provides the SUNDIALS context
     *   (for SUNDIALS 6.0 and newer) to create the solver with.
     */
    SundialsNonlinearSolver
    make_fixed_point_nonlinear_solver(N_Vector y_template,
                                      int      anderson_acceleration_subspace,
                                      const InvocationContext &inv_ctx);
  } // namespace internal
} // namespace SUNDIALS
#endif

DEAL_II_NAMESPACE_CLOSE
#endif
