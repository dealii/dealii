// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2023 by the deal.II authors
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
 * functions in one place.
 *
 * Functions that are not inlined are:
 * - Functions returning PetscErrorCode that are supposed to be called within
 *   PETSc callbacks.
 * - Functions that need access to internal PETSc headers that we don't want
 *   to expose to deal.II users
 */
#ifndef dealii_petsc_compatibility_h
#define dealii_petsc_compatibility_h

#include <deal.II/base/config.h>

#include <deal.II/lac/exceptions.h>

#ifdef DEAL_II_WITH_PETSC

#  include <petscconf.h>
#  include <petscksp.h>
#  include <petscmat.h>
#  include <petscpc.h>
#  include <petscsnes.h>
#  include <petscts.h>

#  include <string>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  /**
   * Set an option in the global PETSc database. This function just wraps
   * PetscOptionsSetValue and checks the error return value.
   */
  inline void
  set_option_value(const std::string &name, const std::string &value)
  {
    const PetscErrorCode ierr =
      PetscOptionsSetValue(nullptr, name.c_str(), value.c_str());
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  /**
   * Set a PETSc matrix option. This function wraps MatSetOption with a
   * version check.
   *
   * @warning The argument option_value is ignored in versions of PETSc
   * before 3.0.0 since the corresponding function did not take this argument.
   */
  inline void
  set_matrix_option(Mat &           matrix,
                    const MatOption option_name,
                    const PetscBool option_value = PETSC_FALSE)
  {
    const PetscErrorCode ierr = MatSetOption(matrix, option_name, option_value);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  /**
   * Tell PETSc that we are not planning on adding new entries to the
   * matrix. Generate errors in debug mode.
   */
  inline void
  close_matrix(Mat &matrix)
  {
#  ifdef DEBUG
    set_matrix_option(matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
#  else
    set_matrix_option(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
#  endif
  }



  /**
   * Tell PETSc to keep the SparsityPattern entries even if we delete a
   * row with clear_rows() which calls MatZeroRows(). Otherwise one can
   * not write into that row afterwards.
   */
  inline void
  set_keep_zero_rows(Mat &matrix)
  {
    set_matrix_option(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  }



  /**
   * Tell PETSc that the status of the vector has changed.
   */
  void
  petsc_increment_state_counter(Vec v);



  /**
   * Tell PETSc that the status of the matrix has changed.
   */
  void
  petsc_increment_state_counter(Mat A);


  /**
   * Set the failed reason for the preconditioner.
   */
  PetscErrorCode
  pc_set_failed_reason(PC pc, PCFailedReason reason);



  /**
   * Resets internal domain error flags in the SNES object.
   */
  void
  snes_reset_domain_flags(SNES snes);



  /**
   * Resets internal domain error flags in the SNES object.
   */
  void
  snes_set_jacobian_domain_error(SNES snes);



  /**
   * Tell PETSc nonlinear solver to use matrix free finite differencing (MFFD).
   *
   * @p mf_operator indicates to use MFFD for the linear system matrix
   * but use a user defined matrix for preconditioning purposed.
   *
   * @p mf indicates to use MFFD for the both the linear system matrix
   * and the preconditioning matrix.
   */
  void
  set_use_matrix_free(SNES snes, const bool mf_operator, const bool mf);



  /**
   * Tell PETSc ODE solver to use matrix free finite differencing (MFFD).
   *
   * @p mf_operator indicates to use MFFD for the linear system matrix
   * but use a user defined matrix for preconditioning purposed.
   *
   * @p mf indicates to use MFFD for the both the linear system matrix
   * and the preconditioning matrix.
   */
  void
  set_use_matrix_free(TS ts, const bool mf_operator, const bool mf);



  /**
   * Set final time for ODE integration.
   */
  void
  ts_set_max_time(TS ts, const PetscReal maxtime);



  /**
   * Set maximum number of steps for ODE integration.
   */
  void
  ts_set_max_steps(TS ts, const PetscInt maxsteps);



  /**
   * Return current step number.
   */
  unsigned int
  ts_get_step_number(TS ts);



  /**
   * Return true if the TS has a SNES object.
   */
  bool
  ts_has_snes(TS ts);

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
#endif // dealii_petsc_compatibility_h
