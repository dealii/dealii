// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_complex_overloads_h
#define dealii_complex_overloads_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
template <typename T>
struct EnableIfScalar;
template <typename T, typename U>
struct ProductType;

#ifndef DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
/**
 * Provide an <tt>operator*</tt> that operates on mixed complex floating point
 * types. Annoyingly, the standard library does not provide such an
 * operator...
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
typename ProductType<std::complex<T>, std::complex<U>>::type inline
operator*(const std::complex<T> &left, const std::complex<U> &right)
{
  typedef
    typename ProductType<std::complex<T>, std::complex<U>>::type result_type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}


/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of a complex
 * floating point type with a different real floating point type. Annoyingly,
 * the standard library does not provide such an operator...
 *
 * @relatesalso EnableIfScalar
 * @relatesalso ProductType
 */
template <typename T, typename U>
typename ProductType<std::complex<T>,
                     typename EnableIfScalar<U>::type>::type inline
operator*(const std::complex<T> &left, const U &right)
{
  typedef typename ProductType<std::complex<T>, U>::type result_type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}


/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of a real
 * floating point type with a different complex floating point type.
 * Annoyingly, the standard library does not provide such an operator...
 *
 * @relatesalso EnableIfScalar
 * @relatesalso ProductType
 */
template <typename T, typename U>
typename ProductType<typename EnableIfScalar<T>::type,
                     std::complex<U>>::type inline
operator*(const T &left, const std::complex<U> &right)
{
  typedef typename ProductType<std::complex<T>, U>::type result_type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}
#endif /* DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS */

DEAL_II_NAMESPACE_CLOSE

#endif
