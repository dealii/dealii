// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

/*
 * Rather than using ifdefs everywhere, try to wrap older versions of PETSc
 * functions in one place.
 */
#ifndef dealii__petsc_compatibility_h
#define dealii__petsc_compatibility_h

#include <deal.II/base/config.h>
#include <deal.II/lac/exceptions.h>

#ifdef DEAL_II_WITH_PETSC

#include <petscconf.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscpc.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
  typedef PetscTruth PetscBooleanType;
#else
  typedef PetscBool PetscBooleanType;
#endif

  /**
   * Set an option in the global PETSc database. This function just wraps
   * PetscOptionsSetValue with a version check (the signature of this function
   * changed in PETSc 3.7.0).
   */
  inline void set_option_value (const std::string &name,
                                const std::string &value)
  {
#if DEAL_II_PETSC_VERSION_LT(3, 7, 0)
    const PetscErrorCode ierr = PetscOptionsSetValue (name.c_str (), value.c_str ());
#else
    const PetscErrorCode ierr = PetscOptionsSetValue (NULL, name.c_str (), value.c_str ());
#endif
    (void)ierr;
    Assert (ierr == 0, ExcPETScError(ierr));
  }



  /**
   * Destroy a PETSc matrix. This function wraps MatDestroy with a version
   * check (the signature of this function changed in PETSc 3.2.0).
   *
   * @warning Since the primary intent of this function is to enable RAII
   * semantics in the PETSc wrappers, this function will not throw an
   * exception if an error occurs, but instead just returns the error code
   * given by MatDestroy.
   *
   */
  inline PetscErrorCode destroy_matrix (Mat &matrix)
  {
    // PETSc will check whether or not matrix is NULL.
#if DEAL_II_PETSC_VERSION_LT(3, 2, 0)
    return MatDestroy (matrix);
#else
    return MatDestroy (&matrix);
#endif
  }



  /**
   * Destroy a Krylov Subspace (KSP) PETSc solver. This function wraps
   * KSPDestroy with a version check (the signature of this function changed
   * in PETSc 3.2.0).
   *
   * @warning Since the primary intent of this function is to enable RAII
   * semantics in the PETSc wrappers, this function will not throw an
   * exception if an error occurs, but instead just returns the error code
   * given by MatDestroy.
   */
  inline PetscErrorCode destroy_krylov_solver (KSP &krylov_solver)
  {
    // PETSc will check whether or not matrix is NULL.
#if DEAL_II_PETSC_VERSION_LT(3, 2, 0)
    return KSPDestroy (krylov_solver);
#else
    return KSPDestroy (&krylov_solver);
#endif
  }



  /**
   * Set a PETSc matrix option. This function wraps MatSetOption with a
   * version check.
   *
   * @warning The argument option_value is ignored in versions of PETSc
   * before 3.0.0 since the corresponding function did not take this argument.
   */
  inline void set_matrix_option (Mat &matrix,
                                 const MatOption option_name,
                                 const PetscBooleanType option_value = PETSC_FALSE)
  {
#if DEAL_II_PETSC_VERSION_LT(3,0,0)
    const int ierr = MatSetOption (matrix, option_name);
#else
    const int ierr = MatSetOption (matrix, option_name, option_value);
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  /**
   * Tell PETSc that we are not planning on adding new entries to the
   * matrix. Generate errors in debug mode.
   */
  inline void close_matrix (Mat &matrix)
  {
#if DEAL_II_PETSC_VERSION_LT(3, 0, 0)
#  ifdef DEBUG
    set_matrix_option (matrix, MAT_NEW_NONZERO_LOCATION_ERR);
#  else
    set_matrix_option (matrix, MAT_NO_NEW_NONZERO_LOCATIONS);
#  endif
#else
#  ifdef DEBUG
    set_matrix_option (matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
#  else
    set_matrix_option (matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
#  endif
#endif
  }



  /**
   * Tell PETSc to keep the SparsityPattern entries even if we delete a
   * row with clear_rows() which calls MatZeroRows(). Otherwise one can
   * not write into that row afterwards.
   */
  inline void set_keep_zero_rows (Mat &matrix)
  {
#if DEAL_II_PETSC_VERSION_LT(3, 1, 0)
    set_matrix_option (matrix, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
#else
    set_matrix_option (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
#endif // dealii__petsc_compatibility_h
