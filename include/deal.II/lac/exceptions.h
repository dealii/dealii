// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
  /** @{ */

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
   * The operation requires a sparsity pattern.
   */
  DeclExceptionMsg(
    ExcNeedsSparsityPattern,
    "This function requires that the current object have a "
    "sparsity pattern attached to it, but no sparsity pattern "
    "is available. This usually means that there is a missing "
    "reinit() call which would have added the sparsity pattern.");

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

  /** @} */
} // namespace LACExceptions


using namespace LACExceptions;


DEAL_II_NAMESPACE_CLOSE

#endif
