// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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

#ifndef dealii__lapack_support_h
#define dealii__lapack_support_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN


namespace LAPACKSupport
{
  /**
   * Most LAPACK functions change the contents of the matrix applied to to
   * something which is not a matrix anymore. Therefore, LAPACK matrix classes
   * in <tt>deal.II</tt> have a state flag indicating what happened to them.
   *
   * @author Guido Kanschat, 2005
   */
  enum State
  {
    /// Contents is actually a matrix.
    matrix,
    /// Contents is the inverse of a matrix.
    inverse_matrix,
    /// Contents is an LU decomposition.
    lu,
    /// Eigenvalue vector is filled
    eigenvalues,
    /// Matrix contains singular value decomposition,
    svd,
    /// Matrix is the inverse of a singular value decomposition
    inverse_svd,
    /// Contents is something useless.
    unusable = 0x8000
  };

  /**
   * Function printing the name of a State.
   */
  inline const char *state_name(State s)
  {
    switch (s)
      {
      case matrix:
        return "matrix";
      case inverse_matrix:
        return "inverse matrix";
      case lu:
        return "lu decomposition";
      case eigenvalues:
        return "eigenvalues";
      case svd:
        return "svd";
      case inverse_svd:
        return "inverse_svd";
      case unusable:
        return "unusable";
      default:
        return "unknown";
      }
  }

  /**
   * A matrix can have certain features allowing for optimization, but hard to
   * test. These are listed here.
   */
  enum Properties
  {
    /// No special properties
    general = 0,
    /// Matrix is symmetric
    symmetric = 1,
    /// Matrix is upper triangular
    upper_triangle = 2,
    /// Matrix is lower triangular
    lower_triangle = 4,
    /// Matrix is diagonal
    diagonal = 6,
    /// Matrix is in upper Hessenberg form
    hessenberg = 8
  };

  /**
   * Character constant.
   */
  static const char A = 'A';
  /**
   * Character constant.
   */
  static const char N = 'N';
  /**
   * Character constant.
   */
  static const char T = 'T';
  /**
   * Character constant.
   */
  static const char U = 'U';
  /**
   * Character constant.
   */
  static const char V = 'V';
  /**
   * Integer constant.
   */
  static const int zero = 0;
  /**
   * Integer constant.
   */
  static const int one = 1;

  /**
   * A LAPACK function returned an error code.
   */
  DeclException2(ExcErrorCode, char *, int,
                 << "The function " << arg1 << " returned with an error code " << arg2);

  /**
   * Exception thrown when a matrix is not in a suitable state for an
   * operation. For instance, a LAPACK routine may have left the matrix in an
   * unusable state, then vmult does not make sense anymore.
   */
  DeclException1(ExcState, State,
                 << "The function cannot be called while the matrix is in state "
                 << state_name(arg1));

  /**
   * This exception is thrown if a certain LAPACK function is not available
   * because no LAPACK installation was detected during configuration.
   */
  DeclException1(ExcMissing, char *,
                 << "When you ran 'cmake' during installation of deal.II, "
                 << "no suitable installation of the BLAS or LAPACK library could "
                 << "be found. Consequently, the function <"
                 << arg1
                 << "> can not be called. Refer to the doc/readme.html "
                 << "file for information on how to ensure that deal.II "
                 << "picks up an existing BLAS and LAPACK installation at "
                 << "configuration time.");
}


DEAL_II_NAMESPACE_CLOSE

#endif
