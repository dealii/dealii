// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

/*
 * Rather than using ifdefs everywhere, try to wrap older versions of PETSc
 * functions in one place. This file contains functions that need access to
 * internal PETSc headers that we don't want to expose to deal.II users
 */
#include <deal.II/lac/petsc_compatibility.h>

#ifdef DEAL_II_WITH_PETSC

#  include <petsc/private/petscimpl.h>
#  include <petsc/private/snesimpl.h>
#  include <petsc/private/tsimpl.h>

// Shorthand notation for PETSc error codes.
#  define AssertPETSc(code)                          \
    do                                               \
      {                                              \
        PetscErrorCode ierr = (code);                \
        AssertThrow(ierr == 0, ExcPETScError(ierr)); \
      }                                              \
    while (0)

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

  void
  snes_reset_domain_flags(SNES snes)
  {
    snes->jacobiandomainerror = PETSC_FALSE;
    snes->domainerror         = PETSC_FALSE;
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
