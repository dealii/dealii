// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__lac_exceptions_h
#define __deal2__lac_exceptions_h

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace LACExceptions
{
  /**
   * @addtogroup Exceptions
   */
  //@{

  /**
   * This function only works for
   * quadratic matrices.
   */
  DeclException0 (ExcNotQuadratic);

  /**
   * The operation cannot be finished
   * since the matrix is singular.
   */
  DeclException0 (ExcSingular);

  /**
   * Block indices of two block
   * objects are different.
   */
  DeclException0 (ExcDifferentBlockIndices);

  /**
   * An error of a PETSc function was
   * encountered. Check the PETSc
   * documentation for details.
   */
  DeclException1 (ExcPETScError,
                  int,
                  << "An error with error number " << arg1
                  << " occurred while calling a PETSc function");

  /**
   * An error of a Trilinos function was
   * encountered. Check the Trilinos
   * documentation for details.
   */
  DeclException1 (ExcTrilinosError,
                  int,
                  << "An error with error number " << arg1
                  << " occurred while calling a Trilinos function");

  //@}
}


using namespace LACExceptions;


DEAL_II_NAMESPACE_CLOSE

#endif
