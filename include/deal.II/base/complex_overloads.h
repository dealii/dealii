// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_complex_overloads_h
#define dealii_complex_overloads_h

#include <deal.II/base/config.h>

#include <complex>
#include <type_traits>


DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN

// Forward declarations
template <typename T, typename U>
struct ProductType;

#endif

#ifndef DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS

/**
 * Provide an <tt>operator*</tt> that operates on mixed complex floating point
 * types. Annoyingly, the standard library does not provide such an
 * operator.
 *
 * @note Because the C++ standard does not provide for mixed-precision complex
 *   operators, code such as the following does not compile:
 *   @code
 *     double              factor = 3.141;
 *     std::complex<float> x(1,2);        // =1+2i
 *     auto z = factor*x;                 // error
 *   @endcode
 *   Because this does not compile, writing mixed-precision complex linear
 * algebra libraries is not easily possible without much additional work that
 * requires explicit casts. For example, one would have to write the code above
 * as follows:
 *   @code
 *     double              factor = 1.0;
 *     std::complex<float> x(1,2);        // =1+2i
 *     using P = ProductType<double,std::complex<float>>::type;
 *     auto z = static_cast<P>(factor) * static_cast<P>(x);
 *   @endcode
 *   This makes much code unreadable. As a consequence, we define the
 *   necessary overloaded multiplication and division operators for mixed
 *   complex arithmetic (in namespace `dealii`) to make the code above
 *   compile without the extra casts.
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
inline std::enable_if_t<
  std::is_floating_point_v<T> && std::is_floating_point_v<U> &&
    !std::is_same_v<T, U>,
  typename ProductType<std::complex<T>, std::complex<U>>::type>
operator*(const std::complex<T> &left, const std::complex<U> &right)
{
  using result_type =
    typename ProductType<std::complex<T>, std::complex<U>>::type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * Provide an <tt>operator/</tt> that operates on mixed complex floating point
 * types. Annoyingly, the standard library does not provide such an
 * operator.
 *
 * @note Because the C++ standard does not provide for mixed-precision complex
 *   operators, code such as the following does not compile:
 *   @code
 *     double              factor = 3.141;
 *     std::complex<float> x(1,2);        // =1+2i
 *     auto z = factor*x;                 // error
 *   @endcode
 *   Because this does not compile, writing mixed-precision complex linear
 * algebra libraries is not easily possible without much additional work that
 * requires explicit casts. For example, one would have to write the code above
 * as follows:
 *   @code
 *     double              factor = 1.0;
 *     std::complex<float> x(1,2);        // =1+2i
 *     using P = ProductType<double,std::complex<float>>::type;
 *     auto z = static_cast<P>(factor) * static_cast<P>(x);
 *   @endcode
 *   This makes much code unreadable. As a consequence, we define the
 *   necessary overloaded multiplication and division operators for mixed
 *   complex arithmetic (in namespace `dealii`) to make the code above
 *   compile without the extra casts.
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
inline std::enable_if_t<
  std::is_floating_point_v<T> && std::is_floating_point_v<U> &&
    !std::is_same_v<T, U>,
  typename ProductType<std::complex<T>, std::complex<U>>::type>
operator/(const std::complex<T> &left, const std::complex<U> &right)
{
  using result_type =
    typename ProductType<std::complex<T>, std::complex<U>>::type;
  return static_cast<result_type>(left) / static_cast<result_type>(right);
}

/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of a complex
 * floating point type with a different real floating point type. Annoyingly,
 * the standard library does not provide such an operator.
 *
 * @note Because the C++ standard does not provide for mixed-precision complex
 *   operators, code such as the following does not compile:
 *   @code
 *     double              factor = 3.141;
 *     std::complex<float> x(1,2);        // =1+2i
 *     auto z = factor*x;                 // error
 *   @endcode
 *   Because this does not compile, writing mixed-precision complex linear
 * algebra libraries is not easily possible without much additional work that
 * requires explicit casts. For example, one would have to write the code above
 * as follows:
 *   @code
 *     double              factor = 1.0;
 *     std::complex<float> x(1,2);        // =1+2i
 *     using P = ProductType<double,std::complex<float>>::type;
 *     auto z = static_cast<P>(factor) * static_cast<P>(x);
 *   @endcode
 *   This makes much code unreadable. As a consequence, we define the
 *   necessary overloaded multiplication and division operators for mixed
 *   complex arithmetic (in namespace `dealii`) to make the code above
 *   compile without the extra casts.
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
inline std::enable_if_t<std::is_floating_point_v<T> &&
                          std::is_floating_point_v<U> && !std::is_same_v<T, U>,
                        typename ProductType<std::complex<T>, U>::type>
operator*(const std::complex<T> &left, const U &right)
{
  using result_type = typename ProductType<std::complex<T>, U>::type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * Provide an <tt>operator/</tt> for a scalar division of a complex
 * floating point type with a different real floating point type. Annoyingly,
 * the standard library does not provide such an operator.
 *
 * @note Because the C++ standard does not provide for mixed-precision complex
 *   operators, code such as the following does not compile:
 *   @code
 *     double              factor = 3.141;
 *     std::complex<float> x(1,2);        // =1+2i
 *     auto z = factor*x;                 // error
 *   @endcode
 *   Because this does not compile, writing mixed-precision complex linear
 * algebra libraries is not easily possible without much additional work that
 * requires explicit casts. For example, one would have to write the code above
 * as follows:
 *   @code
 *     double              factor = 1.0;
 *     std::complex<float> x(1,2);        // =1+2i
 *     using P = ProductType<double,std::complex<float>>::type;
 *     auto z = static_cast<P>(factor) * static_cast<P>(x);
 *   @endcode
 *   This makes much code unreadable. As a consequence, we define the
 *   necessary overloaded multiplication and division operators for mixed
 *   complex arithmetic (in namespace `dealii`) to make the code above
 *   compile without the extra casts.
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
inline std::enable_if_t<std::is_floating_point_v<T> &&
                          std::is_floating_point_v<U> && !std::is_same_v<T, U>,
                        typename ProductType<std::complex<T>, U>::type>
operator/(const std::complex<T> &left, const U &right)
{
  using result_type = typename ProductType<std::complex<T>, U>::type;
  return static_cast<result_type>(left) / static_cast<result_type>(right);
}

/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of a real
 * floating point type with a different complex floating point type.
 * Annoyingly, the standard library does not provide such an operator.
 *
 * @note Because the C++ standard does not provide for mixed-precision complex
 *   operators, code such as the following does not compile:
 *   @code
 *     double              factor = 3.141;
 *     std::complex<float> x(1,2);        // =1+2i
 *     auto z = factor*x;                 // error
 *   @endcode
 *   Because this does not compile, writing mixed-precision complex linear
 * algebra libraries is not easily possible without much additional work that
 * requires explicit casts. For example, one would have to write the code above
 * as follows:
 *   @code
 *     double              factor = 1.0;
 *     std::complex<float> x(1,2);        // =1+2i
 *     using P = ProductType<double,std::complex<float>>::type;
 *     auto z = static_cast<P>(factor) * static_cast<P>(x);
 *   @endcode
 *   This makes much code unreadable. As a consequence, we define the
 *   necessary overloaded multiplication and division operators for mixed
 *   complex arithmetic (in namespace `dealii`) to make the code above
 *   compile without the extra casts.
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
inline std::enable_if_t<std::is_floating_point_v<T> &&
                          std::is_floating_point_v<U> && !std::is_same_v<T, U>,
                        typename ProductType<T, std::complex<U>>::type>
operator*(const T &left, const std::complex<U> &right)
{
  using result_type = typename ProductType<T, std::complex<U>>::type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * Provide an <tt>operator/</tt> for a scalar division of a real
 * floating point type with a different complex floating point type.
 * Annoyingly, the standard library does not provide such an operator.
 *
 * @note Because the C++ standard does not provide for mixed-precision complex
 *   operators, code such as the following does not compile:
 *   @code
 *     double              factor = 3.141;
 *     std::complex<float> x(1,2);        // =1+2i
 *     auto z = factor*x;                 // error
 *   @endcode
 *   Because this does not compile, writing mixed-precision complex linear
 * algebra libraries is not easily possible without much additional work that
 * requires explicit casts. For example, one would have to write the code above
 * as follows:
 *   @code
 *     double              factor = 1.0;
 *     std::complex<float> x(1,2);        // =1+2i
 *     using P = ProductType<double,std::complex<float>>::type;
 *     auto z = static_cast<P>(factor) * static_cast<P>(x);
 *   @endcode
 *   This makes much code unreadable. As a consequence, we define the
 *   necessary overloaded multiplication and division operators for mixed
 *   complex arithmetic (in namespace `dealii`) to make the code above
 *   compile without the extra casts.
 *
 * @relatesalso ProductType
 */
template <typename T, typename U>
inline std::enable_if_t<std::is_floating_point_v<T> &&
                          std::is_floating_point_v<U> && !std::is_same_v<T, U>,
                        typename ProductType<T, std::complex<U>>::type>
operator/(const T &left, const std::complex<U> &right)
{
  using result_type = typename ProductType<T, std::complex<U>>::type;
  return static_cast<result_type>(left) / static_cast<result_type>(right);
}


#endif /* DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS */

DEAL_II_NAMESPACE_CLOSE

#endif
