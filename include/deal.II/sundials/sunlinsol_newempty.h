//-----------------------------------------------------------
//
//    Copyright (C) 2020 - 2021 by the deal.II authors
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

/*
 * The functions in this file are based on an implementation distributed within
 * the SUNDIALS package, see the license here:
 * https://computing.llnl.gov/projects/sundials/license.
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David J. Gardner, Carol S. Woodward, and
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------*/

#ifndef dealii_sundials_sunlinsol_newempty_h
#define dealii_sundials_sunlinsol_newempty_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <sundials/sundials_linearsolver.h>

DEAL_II_NAMESPACE_OPEN
namespace SUNDIALS
{
  namespace internal
  {
    /**
     * Create a new SUNLinearSolver structure without any content and
     * operations set to `nullptr`.
     */
    inline SUNLinearSolver
    SUNLinSolNewEmpty()
    {
      /* create linear solver object */
      SUNLinearSolver LS = new _generic_SUNLinearSolver;

      /* create linear solver ops structure */
      SUNLinearSolver_Ops ops = new _generic_SUNLinearSolver_Ops;

      /* initialize operations to nullptr */
      ops->gettype           = nullptr;
      ops->setatimes         = nullptr;
      ops->setpreconditioner = nullptr;
      ops->setscalingvectors = nullptr;
      ops->initialize        = nullptr;
      ops->setup             = nullptr;
      ops->solve             = nullptr;
      ops->numiters          = nullptr;
      ops->resnorm           = nullptr;
      ops->resid             = nullptr;
      ops->lastflag          = nullptr;
      ops->space             = nullptr;
      ops->free              = nullptr;

      /* attach ops and initialize content to nullptr */
      LS->ops     = ops;
      LS->content = nullptr;

      return (LS);
    }

    /**
     * Free the memory associated with @p solver which was previously allocated
     * with a call to SUNLinSolNewEmpty().
     *
     * @note A call to this function does not deallocate the `content` field.
     *
     * @param solver The solver memory to free
     */
    inline void
    SUNLinSolFreeEmpty(SUNLinearSolver solver)
    {
      if (solver == nullptr)
        return;

      /* free non-nullptr ops structure */
      if (solver->ops)
        delete solver->ops;
      solver->ops = nullptr;

      /* free overall linear solver object */
      delete solver;
      return;
    }

  } // namespace internal
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS
#endif // dealii_sundials_sunlinsol_newempty_h
