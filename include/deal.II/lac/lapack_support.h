// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_lapack_support_h
#define dealii_lapack_support_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>

DEAL_II_NAMESPACE_OPEN

namespace types
{
#ifdef LAPACK_WITH_64BIT_BLAS_INDICES
  /**
   * Integer type in BLAS.
   */
  using blas_int = long long;
#else
  /**
   * Integer type in BLAS.
   */
  using blas_int = int;
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
   * %Function printing the name of a State.
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
   * %Function printing the name of a Property.
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

    DEAL_II_NOT_IMPLEMENTED();
    return "invalid";
  }

  /**
   * Character constant.
   */
  constexpr char A = 'A';
  /**
   * Character constant.
   */
  constexpr char N = 'N';
  /**
   * Character constant.
   */
  constexpr char O = 'O';
  /**
   * Character constant.
   */
  constexpr char T = 'T';
  /**
   * Character constant.
   */
  constexpr char U = 'U';
  /**
   * Character constant.
   */
  constexpr char L = 'L';
  /**
   * Character constant.
   */
  constexpr char V = 'V';
  /**
   * Integer constant.
   */
  constexpr types::blas_int zero = 0;
  /**
   * Integer constant.
   */
  constexpr types::blas_int one = 1;

  /**
   * A LAPACK function returned an error code.
   */
  DeclException2(ExcErrorCode,
                 std::string,
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
    std::string,
    << "When you ran 'cmake' during installation of deal.II, no suitable "
    << "installation of the BLAS or LAPACK library could be found. "
    << "Consequently, the function <" << arg1 << "> can not be called. "
    << "Refer to the readme at https://dealii.org/current/readme.html for "
    << "information on how to ensure that deal.II picks up an existing "
    << "BLAS and LAPACK installation at configuration time.");
} // namespace LAPACKSupport


DEAL_II_NAMESPACE_CLOSE

#endif
