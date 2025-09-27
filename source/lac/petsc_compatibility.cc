// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Rather than using ifdefs everywhere, try to wrap older versions of PETSc
 * functions in one place. This file contains functions that need access to
 * internal PETSc headers that we don't want to expose to deal.II users
 */
#include <deal.II/lac/petsc_compatibility.h>

#ifdef DEAL_II_WITH_PETSC

#  include <petsc/private/pcimpl.h>
#  include <petsc/private/petscimpl.h>
#  include <petsc/private/snesimpl.h>
#  include <petsc/private/tsimpl.h>
#  include <petscdm.h>

// Shorthand notation for PETSc error codes.
#  define AssertPETSc(code)                          \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (false)

DEAL_II_NAMESPACE_OPEN


namespace PETScWrappers
{
  void
  petsc_increment_state_counter(Vec v)
  {
    AssertPETSc(PetscObjectStateIncrease(reinterpret_cast<PetscObject>(v)));
  }

  void
  petsc_increment_state_counter(Mat A)
  {
    AssertPETSc(PetscObjectStateIncrease(reinterpret_cast<PetscObject>(A)));
  }

  PetscErrorCode
  pc_set_failed_reason(PC pc, PCFailedReason reason)
  {
#  if DEAL_II_PETSC_VERSION_GTE(3, 14, 0)
    return PCSetFailedReason(pc, reason);
#  else
    pc->failedreason = reason;
    return 0;
#  endif
  }

  void
  snes_reset_domain_flags(SNES snes)
  {
#  if DEAL_II_PETSC_VERSION_GTE(3, 11, 0)
    snes->jacobiandomainerror = PETSC_FALSE;
#  endif
    snes->domainerror = PETSC_FALSE;
  }

  void
  snes_set_jacobian_domain_error(SNES snes)
  {
#  if DEAL_II_PETSC_VERSION_GTE(3, 11, 0)
    snes->jacobiandomainerror = PETSC_TRUE;
#  else
    // There is no equivalent, and since this used to stop
    // computations, we opt to set the converged reason
    snes->reason = SNES_DIVERGED_FUNCTION_DOMAIN;
#  endif
  }

  void
  set_use_matrix_free(SNES snes, bool mf_operator, bool mf)
  {
#  if DEAL_II_PETSC_VERSION_LT(3, 13, 1)
    snes->mf          = mf ? PETSC_TRUE : PETSC_FALSE;
    snes->mf_operator = mf_operator ? PETSC_TRUE : PETSC_FALSE;
#  else
    AssertPETSc(SNESSetUseMatrixFree(snes,
                                     mf_operator ? PETSC_TRUE : PETSC_FALSE,
                                     mf ? PETSC_TRUE : PETSC_FALSE));
#  endif
  }

  void
  set_use_matrix_free(TS ts, const bool mf_operator, const bool mf)
  {
    SNES snes;
    AssertPETSc(TSGetSNES(ts, &snes));
    set_use_matrix_free(snes, mf_operator, mf);
  }

  void
  ts_set_max_steps(TS ts, const PetscInt maxsteps)
  {
#  if DEAL_II_PETSC_VERSION_LT(3, 8, 0)
    if (maxsteps >= 0)
      ts->max_steps = maxsteps;
#  else
    AssertPETSc(TSSetMaxSteps(ts, maxsteps));
#  endif
  }

  void
  ts_set_max_time(TS ts, const PetscReal maxtime)
  {
#  if DEAL_II_PETSC_VERSION_LT(3, 8, 0)
    if (maxtime != PETSC_DEFAULT)
      ts->max_time = maxtime;
#  else
    AssertPETSc(TSSetMaxTime(ts, maxtime));
#  endif
  }

  void
  ts_reset_dm(TS ts)
  {
    AssertPETSc(DMDestroy(&ts->dm));
  }

  unsigned int
  ts_get_step_number(TS ts)
  {
    PetscInt step;
#  if DEAL_II_PETSC_VERSION_LT(3, 8, 0)
    step = ts->steps;
#  else
    AssertPETSc(TSGetStepNumber(ts, &step));
#  endif
    return static_cast<unsigned int>(step);
  }

  bool
  ts_has_snes(TS ts)
  {
    return ts->snes ? true : false;
  }

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
