// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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


// test mixed precision product types between vectorized arrays and floating
// point types

#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <typename Number1,
          std::size_t width1,
          typename Number2,
          std::size_t width2>
using VPT = typename internal::
  VectorizationProductTypeImpl<Number1, width1, Number2, width2>;


template <typename Number1,
          std::size_t width1,
          typename Number2,
          std::size_t width2>
using product_t =
  VectorizedArray<typename VPT<Number1, width1, Number2, width2>::value_type,
                  VPT<Number1, width1, Number2, width2>::width>;


// Test 1: LHS = VectorizedArray; RHS = Scalar type
template <std::size_t width_float, std::size_t width_double>
void
do_test_vectorized_array_and_scalar()
{
  static_assert(width_float >= width_double, "Invalid vectorization width");

  using float_t             = float;
  using double_t            = double;
  using complex_float_t     = std::complex<float_t>;
  using complex_double_t    = std::complex<double_t>;
  using va_float_t          = VectorizedArray<float_t, width_float>;
  using va_double_t         = VectorizedArray<double_t, width_double>;
  using va_complex_float_t  = VectorizedArray<complex_float_t, width_float>;
  using va_complex_double_t = VectorizedArray<complex_double_t, width_double>;

  constexpr const std::size_t width_scalar = 1;

  // Exactly the same types in the vectorized array and the other type
  static_assert(
    std::is_same<product_t<float_t, width_float, float_t, width_scalar>,
                 va_float_t>::value,
    "Not same");
  static_assert(
    std::is_same<product_t<double_t, width_double, double_t, width_scalar>,
                 va_double_t>::value,
    "Not same");

  static_assert(
    std::is_same<
      product_t<complex_float_t, width_float, complex_float_t, width_scalar>,
      va_complex_float_t>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<complex_double_t, width_double, complex_double_t, width_scalar>,
      va_complex_double_t>::value,
    "Not same");

  // Different types in the vectorized array and the other type
  // - Not mixing complex and non-complex types
  static_assert(
    std::is_same<product_t<double_t, width_double, float_t, width_scalar>,
                 va_double_t>::value,
    "Not same");
  static_assert(
    std::is_same<product_t<float_t, width_float, double_t, width_scalar>,
                 typename std::conditional<width_float == width_double,
                                           va_double_t,
                                           va_float_t>::type>::value,
    "Not same");

  static_assert(
    std::is_same<
      product_t<complex_double_t, width_double, complex_float_t, width_scalar>,
      va_complex_double_t>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<complex_float_t, width_float, complex_double_t, width_scalar>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");

  // Different types in the vectorized array and the other type
  // - Mixing complex and non-complex types
  static_assert(
    std::is_same<
      product_t<complex_double_t, width_double, float_t, width_scalar>,
      va_complex_double_t>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<double_t, width_double, complex_float_t, width_scalar>,
      va_complex_double_t>::value,
    "Not same");

  static_assert(
    std::is_same<
      product_t<complex_float_t, width_float, double_t, width_scalar>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<float_t, width_float, complex_double_t, width_scalar>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");
}


// Test 2: LHS = VectorizedArray; RHS = VectorizedArray
template <std::size_t width_float, std::size_t width_double>
void
do_test_vectorized_array_and_vectorized_array()
{
  static_assert(width_float >= width_double, "Invalid vectorization width");

  using float_t             = float;
  using double_t            = double;
  using complex_float_t     = std::complex<float_t>;
  using complex_double_t    = std::complex<double_t>;
  using va_float_t          = VectorizedArray<float_t, width_float>;
  using va_double_t         = VectorizedArray<double_t, width_double>;
  using va_complex_float_t  = VectorizedArray<complex_float_t, width_float>;
  using va_complex_double_t = VectorizedArray<complex_double_t, width_double>;

  // Exactly the same types in the vectorized array and the other type
  static_assert(
    std::is_same<product_t<float_t, width_float, float_t, width_float>,
                 va_float_t>::value,
    "Not same");
  static_assert(
    std::is_same<product_t<double_t, width_double, double_t, width_double>,
                 va_double_t>::value,
    "Not same");

  static_assert(
    std::is_same<
      product_t<complex_float_t, width_float, complex_float_t, width_float>,
      va_complex_float_t>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<complex_double_t, width_double, complex_double_t, width_double>,
      va_complex_double_t>::value,
    "Not same");

  // Different types in the vectorized array and the other type
  // - Not mixing complex and non-complex types
  static_assert(
    std::is_same<product_t<double_t, width_double, float_t, width_float>,
                 typename std::conditional<width_float == width_double,
                                           va_double_t,
                                           va_float_t>::type>::value,
    "Not same");
  static_assert(
    std::is_same<product_t<float_t, width_float, double_t, width_double>,
                 typename std::conditional<width_float == width_double,
                                           va_double_t,
                                           va_float_t>::type>::value,
    "Not same");

  static_assert(
    std::is_same<
      product_t<complex_double_t, width_double, complex_float_t, width_float>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<complex_float_t, width_float, complex_double_t, width_double>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");

  // Different types in the vectorized array and the other type
  // - Mixing complex and non-complex types
  static_assert(
    std::is_same<
      product_t<complex_double_t, width_double, float_t, width_float>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<double_t, width_double, complex_float_t, width_float>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");

  static_assert(
    std::is_same<
      product_t<complex_float_t, width_float, double_t, width_double>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");
  static_assert(
    std::is_same<
      product_t<float_t, width_float, complex_double_t, width_double>,
      typename std::conditional<width_float == width_double,
                                va_complex_double_t,
                                va_complex_float_t>::type>::value,
    "Not same");
}


int
main()
{
  initlog();

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test_vectorized_array_and_scalar<16, 8>();
  do_test_vectorized_array_and_vectorized_array<16, 8>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test_vectorized_array_and_scalar<8, 4>();
  do_test_vectorized_array_and_vectorized_array<8, 4>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test_vectorized_array_and_scalar<4, 2>();
  do_test_vectorized_array_and_vectorized_array<4, 2>();
#endif

  do_test_vectorized_array_and_scalar<1, 1>();
  do_test_vectorized_array_and_vectorized_array<1, 1>();

  deallog << "OK" << std::endl;
}
