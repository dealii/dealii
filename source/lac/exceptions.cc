// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/lac/exceptions.h>

#ifdef DEAL_II_WITH_PETSC
#  include <petscconf.h>
#  include <petscsys.h>
#endif // DEAL_II_WITH_PETSC

DEAL_II_NAMESPACE_OPEN

namespace LACExceptions
{
  ExcPETScError::ExcPETScError(const int error_code) : error_code(error_code)
  {}

  void
  ExcPETScError::print_info(std::ostream &out) const
  {
    out << "deal.II encountered an error while calling a PETSc function."
        << std::endl;
#ifdef DEAL_II_WITH_PETSC
    // PetscErrorMessage changes the value in a pointer to refer to a
    // statically allocated description of the current error message.
    const char *         petsc_message;
    const PetscErrorCode ierr = PetscErrorMessage(error_code,
                                                  &petsc_message,
                                                  /*specific=*/nullptr);
    if (ierr == 0 && petsc_message != nullptr)
      {
        out << "The description of the error provided by PETSc is \""
            << petsc_message << "\"." << std::endl;
      }
    else
      {
        out
          << "PETSc was not able to determine a description for this particular error code."
          << std::endl;
      }
#endif // DEAL_II_WITH_PETSC
    out << "The numerical value of the original error code is " << error_code
        << "." << std::endl;
  }
} // namespace LACExceptions

DEAL_II_NAMESPACE_CLOSE
