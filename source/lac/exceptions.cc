// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/lac/exceptions.h>

#include <iostream>

#ifdef DEAL_II_WITH_PETSC
#  include <petscsys.h>
#endif // DEAL_II_WITH_PETSC


DEAL_II_NAMESPACE_OPEN

namespace LACExceptions
{
  ExcPETScError::ExcPETScError(const int error_code)
    : error_code(error_code)
  {}

  void
  ExcPETScError::print_info(std::ostream &out) const
  {
    out << "deal.II encountered an error while calling a PETSc function."
        << std::endl;
#ifdef DEAL_II_WITH_PETSC
    // PetscErrorMessage changes the value in a pointer to refer to a
    // statically allocated description of the current error message.
    const char          *petsc_message;
    const PetscErrorCode ierr =
      PetscErrorMessage(static_cast<PetscErrorCode>(error_code),
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
        << '.' << std::endl;
  }
} // namespace LACExceptions

DEAL_II_NAMESPACE_CLOSE
