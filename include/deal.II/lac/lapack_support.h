// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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

#ifndef dealii_lapack_support_h
#define dealii_lapack_support_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace types
{
#ifdef LAPACK_WITH_64BIT_BLAS_INDICES
  /**
   * Integer type in BLAS.
   */
  typedef long long blas_int;
#else
  /**
   * Integer type in BLAS.
   */
  typedef int blas_int;
#endif
} // namespace types

/**
 * A namespace containing constants, exceptions, enumerations, and other
 * utilities used by the deal.II LAPACK bindings.
 */
namespace LAPACKSupport
{
  /**
   * Most of the LAPACK functions one can apply to a matrix (e.g., by calling
   * the member functions of this class) change its content in some ways. For
   * example, they may invert the matrix, or may replace it by a matrix whose
   * columns represent the eigenvectors of the original content of the matrix.
   * The elements of this enumeration are therefore used to track what is
   * currently being stored by this object.
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
    /// Contents is a Cholesky decomposition.
    cholesky,
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
  inline const char *
  state_name(State s)
  {
    switch (s)
      {
        case matrix:
          return "matrix";
        case inverse_matrix:
          return "inverse matrix";
        case lu:
          return "lu decomposition";
        case cholesky:
          return "cholesky decomposition";
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
  enum Property
  {
    /// No special properties
    general = 0,
    /// Matrix is symmetric
    symmetric = 1,
    /// Matrix is upper triangular
    upper_triangular = 2,
    /// Matrix is lower triangular
    lower_triangular = 4,
    /// Matrix is diagonal
    diagonal = 6,
    /// Matrix is in upper Hessenberg form
    hessenberg = 8
  };

  /**
   * Function printing the name of a Property.
   */
  inline const char *
  property_name(const Property s)
  {
    switch (s)
      {
        case general:
          return "general";
        case symmetric:
          return "symmetric";
        case upper_triangular:
          return "upper triangular";
        case lower_triangular:
          return "lower triangular";
        case diagonal:
          return "diagonal";
        case hessenberg:
          return "Hessenberg";
      }

    Assert(false, ExcNotImplemented());
    return "invalid";
  }

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
  static const char L = 'L';
  /**
   * Character constant.
   */
  static const char V = 'V';
  /**
   * Integer constant.
   */
  static const types::blas_int zero = 0;
  /**
   * Integer constant.
   */
  static const types::blas_int one = 1;

  /**
   * A LAPACK function returned an error code.
   */
  DeclException2(ExcErrorCode,
                 char *,
                 types::blas_int,
                 << "The function " << arg1 << " returned with an error code "
                 << arg2);

  /**
   * Exception thrown when a matrix is not in a suitable state for an
   * operation. For instance, a LAPACK routine may have left the matrix in an
   * unusable state, then vmult does not make sense anymore.
   */
  DeclException1(
    ExcState,
    State,
    << "The function cannot be called while the matrix is in state "
    << state_name(arg1));

  /**
   * Exception thrown when a matrix does not have suitable properties for an
   * operation.
   */
  DeclException1(ExcProperty,
                 Property,
                 << "The function cannot be called with a "
                 << property_name(arg1) << " matrix.");

  /**
   * This exception is thrown if a certain LAPACK function is not available
   * because no LAPACK installation was detected during configuration.
   */
  DeclException1(
    ExcMissing,
    char *,
    << "When you ran 'cmake' during installation of deal.II, "
    << "no suitable installation of the BLAS or LAPACK library could "
    << "be found. Consequently, the function <" << arg1
    << "> can not be called. Refer to the doc/readme.html "
    << "file for information on how to ensure that deal.II "
    << "picks up an existing BLAS and LAPACK installation at "
    << "configuration time.");
} // namespace LAPACKSupport


DEAL_II_NAMESPACE_CLOSE

#endif
