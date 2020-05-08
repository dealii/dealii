// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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

#ifndef dealii_lac_exceptions_h
#define dealii_lac_exceptions_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace LACExceptions
{
  /**
   * @addtogroup Exceptions
   */
  //@{

  /**
   * This function only works for quadratic matrices.
   */
  DeclExceptionMsg(ExcNotQuadratic,
                   "This function only works for quadratic objects!");

  /**
   * The operation cannot be finished since the matrix is singular.
   */
  DeclException0(ExcSingular);

  /**
   * Block indices of two block objects are different.
   */
  DeclException0(ExcDifferentBlockIndices);

  /**
   * Exception thrown when a PETSc function reports an error. If possible,
   * this exception uses the message provided by
   * <code>PetscErrorMessage</code> to print a description of the error.
   *
   * @note For backwards compatibility this is defined whether or not deal.II
   * is compiled with PETSc.
   */
  class ExcPETScError : public dealii::ExceptionBase
  {
  public:
    ExcPETScError(const int error_code);

    virtual void
    print_info(std::ostream &out) const override;

    const int error_code;
  };

  /**
   * An error of a Trilinos function was encountered. Check the Trilinos
   * documentation for details.
   */
  DeclException1(ExcTrilinosError,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while calling a Trilinos function");

  //@}
} // namespace LACExceptions


using namespace LACExceptions;


DEAL_II_NAMESPACE_CLOSE

#endif
