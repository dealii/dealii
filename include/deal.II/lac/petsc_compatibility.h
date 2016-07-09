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

#ifdef DEAL_II_WITH_PETSC

#include <petscconf.h>
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
    PetscOptionsSetValue (name.c_str (), value.c_str ());
#else
    PetscOptionsSetValue (NULL, name.c_str (), value.c_str ());
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
#endif // dealii__petsc_compatibility_h
