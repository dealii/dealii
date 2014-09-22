// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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


#ifndef __deal2__vectorization_h
#define __deal2__vectorization_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <cmath>

// Note:
// The flag DEAL_II_COMPILER_VECTORIZATION_LEVEL is essentially constructed
// according to the following scheme
// #ifdef __AVX512F__
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 3
// #ifdef __AVX__
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 2
// #elif defined (__SSE2__)
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 1
// #else
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 0
// #endif
// In addition to checking the flags __AVX__ and __SSE2__, a configure test
// ensures that these feature are not only present but also working properly.

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 // AVX, AVX-512
#include <immintrin.h>
#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL == 1 // SSE2
#include <emmintrin.h>
#endif



// forward declarations
namespace dealii
{
  template <typename Number> class VectorizedArray;
}
namespace std
{
  template <typename Number> ::dealii::VectorizedArray<Number>
  sqrt(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  abs(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  max(const ::dealii::VectorizedArray<Number> &, const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  min (const ::dealii::VectorizedArray<Number> &, const ::dealii::VectorizedArray<Number> &);
}


DEAL_II_NAMESPACE_OPEN


// for safety, also check that __AVX512F__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL == 3  && defined(__AVX512F__)

/**
 * Specialization of VectorizedArray class for double and AVX-512.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 8;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const double x)
  {
    data = _mm512_set_pd(x, x, x, x, x, x, x, x);
    return *this;
  }

  /**
   * Access operator.
   */
  double &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 8);
    return *(reinterpret_cast<double *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 8);
    return *(reinterpret_cast<const double *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm512_add_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm512_sub_pd(data,vec.data);
#endif
    return *this;
  }
  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm512_mul_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm512_div_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 64 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  void load (const double *ptr)
  {
    data = _mm512_loadu_pd (ptr);
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 64 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  void store (double *ptr) const
  {
    _mm512_storeu_pd (ptr, data);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m512d data;

private:
  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm512_sqrt_pd(data);
    return res;
  }

  /**
   * Returns the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs () const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m256d mask = _mm512_set_pd (-0., -0., -0., -0., -0., -0., -0., -0.);
    VectorizedArray res;
    res.data = _mm512_andnot_pd(mask, data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_max_pd (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_pd (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};



/**
 * Specialization for float and AVX.
 */
template<>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 16;

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator = (const float x)
  {
    data = _mm512_set_ps(x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x);
    return *this;
  }

  /**
   * Access operator.
   */
  float &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 16);
    return *(reinterpret_cast<float *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const float &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 16);
    return *(reinterpret_cast<const float *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm512_add_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm512_sub_ps(data,vec.data);
#endif
    return *this;
  }
  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm512_mul_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm512_div_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 64 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  void load (const float *ptr)
  {
    data = _mm512_loadu_ps (ptr);
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 64 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  void store (float *ptr) const
  {
    _mm512_storeu_ps (ptr, data);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m512 data;

private:

  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm512_sqrt_ps(data);
    return res;
  }

  /**
   * Returns the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs () const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m256 mask = _mm512_setzero_ps (-0.f, -0.f, -0.f, -0.f, -0.f, -0.f, -0.f, -0.f,
                                     -0.f, -0.f, -0.f, -0.f, -0.f, -0.f, -0.f, -0.f);
    VectorizedArray res;
    res.data = _mm512_andnot_ps(mask, data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_max_ps (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_ps (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};



#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL == 2  && defined(__AVX__)

/**
 * Specialization of VectorizedArray class for double and AVX.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected
   * in this class.
   */
  static const unsigned int n_array_elements = 4;

  /**
   * This function can be used to set all data
   * fields to a given scalar.
   */
  VectorizedArray &
  operator = (const double x)
  {
    data = _mm256_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  double &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 4);
    return *(reinterpret_cast<double *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 4);
    return *(reinterpret_cast<const double *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
    // if the compiler supports vector
    // arithmetics, we can simply use += operator
    // on the given data type. Otherwise, we need
    // to use the built-in intrinsic command for
    // __m256d
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm256_add_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm256_sub_pd(data,vec.data);
#endif
    return *this;
  }
  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm256_mul_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm256_div_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 32 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  void load (const double *ptr)
  {
    data = _mm256_loadu_pd (ptr);
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 32 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  void store (double *ptr) const
  {
    _mm256_storeu_pd (ptr, data);
  }

  /**
   * Actual data field. Since this class
   * represents a POD data type, it remains
   * public.
   */
  __m256d data;

private:
  /**
   * Returns the square root of this field. Not
   * for use in user code. Use sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_pd(data);
    return res;
  }

  /**
   * Returns the absolute value of this
   * field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs () const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m256d mask = _mm256_set1_pd (-0.);
    VectorizedArray res;
    res.data = _mm256_andnot_pd(mask, data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this
   * field and another one. Not for use in user
   * code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_max_pd (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this
   * field and another one. Not for use in user
   * code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_pd (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};



/**
 * Specialization for float and AVX.
 */
template<>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected
   * in this class.
   */
  static const unsigned int n_array_elements = 8;

  /**
   * This function can be used to set all data
   * fields to a given scalar.
   */
  VectorizedArray &
  operator = (const float x)
  {
    data = _mm256_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  float &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 8);
    return *(reinterpret_cast<float *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const float &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 8);
    return *(reinterpret_cast<const float *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm256_add_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm256_sub_ps(data,vec.data);
#endif
    return *this;
  }
  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm256_mul_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm256_div_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 32 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  void load (const float *ptr)
  {
    data = _mm256_loadu_ps (ptr);
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 32 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  void store (float *ptr) const
  {
    _mm256_storeu_ps (ptr, data);
  }

  /**
   * Actual data field. Since this class
   * represents a POD data type, it remains
   * public.
   */
  __m256 data;

private:

  /**
   * Returns the square root of this field. Not
   * for use in user code. Use sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_ps(data);
    return res;
  }
  /**
   * Returns the absolute value of this
   * field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs () const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m256 mask = _mm256_set1_ps (-0.f);
    VectorizedArray res;
    res.data = _mm256_andnot_ps(mask, data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this
   * field and another one. Not for use in user
   * code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_max_ps (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this
   * field and another one. Not for use in user
   * code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_ps (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};



// for safety, also check that __SSE2__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 1 && defined(__SSE2__)



/**
 * Specialization for double and SSE2.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected
   * in this class.
   */
  static const unsigned int n_array_elements = 2;

  /**
   * This function can be used to set all data
   * fields to a given scalar.
   */

  VectorizedArray &
  operator = (const double x)
  {
    data = _mm_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  double &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 2);
    return *(reinterpret_cast<double *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 2);
    return *(reinterpret_cast<const double *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm_add_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm_sub_pd(data,vec.data);
#endif
    return *this;
  }
  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm_mul_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm_div_pd(data,vec.data);
#endif
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  void load (const double *ptr)
  {
    data = _mm_loadu_pd (ptr);
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  void store (double *ptr) const
  {
    _mm_storeu_pd (ptr, data);
  }

  /**
   * Actual data field. Since this class
   * represents a POD data type, it remains
   * public.
   */
  __m128d data;

private:
  /**
   * Returns the square root of this field. Not
   * for use in user code. Use sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_pd(data);
    return res;
  }

  /**
   * Returns the absolute value of this
   * field. Not for use in user code. Use abs(x)
   * instead.
   */
  VectorizedArray
  get_abs () const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m128d mask = _mm_set1_pd (-0.);
    VectorizedArray res;
    res.data = _mm_andnot_pd(mask, data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this
   * field and another one. Not for use in user
   * code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_max_pd (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this
   * field and another one. Not for use in user
   * code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_min_pd (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};



/**
 * Specialization for float and SSE2.
 */
template <>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected
   * in this class.
   */
  static const unsigned int n_array_elements = 4;

  /**
   * This function can be used to set all data
   * fields to a given scalar.
   */

  VectorizedArray &
  operator = (const float x)
  {
    data = _mm_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  float &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 4);
    return *(reinterpret_cast<float *>(&data)+comp);
  }

  /**
   * Constant access operator.
   */
  const float &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 4);
    return *(reinterpret_cast<const float *>(&data)+comp);
  }

  /**
   * Addition.
   */
  VectorizedArray &
  operator += (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#else
    data = _mm_add_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Subtraction.
   */
  VectorizedArray &
  operator -= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#else
    data = _mm_sub_ps(data,vec.data);
#endif
    return *this;
  }
  /**
   * Multiplication.
   */
  VectorizedArray &
  operator *= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#else
    data = _mm_mul_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Division.
   */
  VectorizedArray &
  operator /= (const VectorizedArray &vec)
  {
#ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#else
    data = _mm_div_ps(data,vec.data);
#endif
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  void load (const float *ptr)
  {
    data = _mm_loadu_ps (ptr);
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  void store (float *ptr) const
  {
    _mm_storeu_ps (ptr, data);
  }

  /**
   * Actual data field. Since this class
   * represents a POD data type, it remains
   * public.
   */
  __m128 data;

private:
  /**
   * Returns the square root of this field. Not
   * for use in user code. Use sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_ps(data);
    return res;
  }

  /**
   * Returns the absolute value of this
   * field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs () const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m128 mask = _mm_set1_ps (-0.f);
    VectorizedArray res;
    res.data = _mm_andnot_ps(mask, data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this
   * field and another one. Not for use in user
   * code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_max_ps (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this
   * field and another one. Not for use in user
   * code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_min_ps (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};


#endif // if DEAL_II_COMPILER_VECTORIZATION_LEVEL > 0


/**
 * This generic class defines a unified interface to a vectorized data
 * type. For general template arguments, this class simply corresponds to
 * the template argument. For example, VectorizedArray<long double> is
 * nothing else but a wrapper around <tt>long double</tt> with exactly one
 * data field of type <tt>long double</tt> and overloaded arithmetic
 * operations. This means that <tt>VectorizedArray<ComplicatedType></tt> has
 * a similar layout as ComplicatedType, provided that ComplicatedType
 * defines basic arithmetic operations. For floats and doubles, an array of
 * numbers are packed together, though. The number of elements packed
 * together depend on the computer system and compiler flags that are used
 * for compilation of deal.II. The fundamental idea of these packed data
 * types is to use one single CPU instruction to perform arithmetic
 * operations on the whole array using the processor's vector units. Most
 * computer systems by 2010 standards will use an array of two doubles and
 * four floats, respectively (this corresponds to the SSE/SSE2 data sets)
 * when compiling deal.II on 64-bit operating systems. On Intel Sandy Bridge
 * processors and newer or AMD Bulldozer processors and newer, four doubles
 * and eight floats are used when deal.II is configured e.g. using gcc with
 * --with-cpu=native or --with-cpu=corei7-avx.
 *
 * This behavior of this class is made similar to the basic data types
 * double and float. The definition of a vectorized array does not
 * initialize the data field but rather leaves it undefined, as is the case
 * for double and float. However, when calling something like
 * VectorizedArray<double> a = VectorizedArray<double>(), it sets all numbers in this
 * field to zero. In other words, this class is a plain old data (POD) type
 * which has an equivalent C representation and can e.g. be safely copied
 * with std::memcpy. This POD layout is also necessary for ensuring correct
 * alignment of data with address boundaries when collected in a vector
 * (i.e., when the first element in a vector is properly aligned, all
 * subsequent elements will be correctly aligned, too).
 *
 * Note that for proper functioning of this class, certain data alignment
 * rules must be respected. This is because the computer expects the
 * starting address of a VectorizedArray<double> field at specific addresses
 * in memory (usually, the address of the vectorized array should be a
 * multiple of the length of the array in bytes). Otherwise, a segmentation
 * fault or a severe loss of performance might occur. When creating a single
 * data field on the stack like <tt>VectorizedArray<double> a =
 * VectorizedArray<double>()</tt>, the compiler will take care of data
 * alignment automatically. However, when allocating a long vector of
 * VectorizedArray<double> data, one needs to respect these rules. Use the
 * class AlignedVector for this purpose. It is a class very similar to
 * std::vector otherwise but always makes sure that data is correctly
 * aligned.
 *
 * @author Katharina Kormann, Martin Kronbichler, 2010, 2011
 */
template <typename Number>
class VectorizedArray
{
public:
  /**
   * This gives the number of vectors collected
   * in this class.
   */
  static const unsigned int n_array_elements = 1;

  // POD means that there should be no
  // user-defined constructors, destructors and
  // copy functions (the standard is somewhat
  // relaxed in C++2011, though).

  /**
   * This function assigns a scalar to this
   * class.
   */

  VectorizedArray &
  operator = (const Number scalar)
  {
    data = scalar;
    return *this;
  }

  /**
   * Access operator (only valid with component
   * 0)
   */
  Number &
  operator [] (const unsigned int comp)
  {
    AssertIndexRange (comp, 1);
    return data;
  }

  /**
   * Constant access operator (only valid with
   * component 0)
   */
  const Number &
  operator [] (const unsigned int comp) const
  {
    AssertIndexRange (comp, 1);
    return data;
  }

  /**
   * Addition
   */
  VectorizedArray &
  operator += (const VectorizedArray<Number> &vec)
  {
    data+=vec.data;
    return *this;
  }

  /**
   * Subtraction
   */
  VectorizedArray &
  operator -= (const VectorizedArray<Number> &vec)
  {
    data-=vec.data;
    return *this;
  }

  /**
   * Multiplication
   */
  VectorizedArray &
  operator *= (const VectorizedArray<Number> &vec)
  {
    data*=vec.data;
    return *this;
  }

  /**
   * Division
   */
  VectorizedArray &
  operator /= (const VectorizedArray<Number> &vec)
  {
    data/=vec.data;
    return *this;
  }

  /**
   * Loads @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by the amount of bytes
   * in the vectorized array, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  void load (const Number *ptr)
  {
    data = *ptr;
  }

  /**
   * Writes the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * the amount of bytes in the vectorized array, as opposed to casting a
   * double address to VectorizedArray<double>*.
   */
  void store (Number *ptr) const
  {
    *ptr = data;
  }

  /**
   * Actual data field. Since this class
   * represents a POD data type, it is declared
   * public.
   */
  Number data;

private:
  /**
   * Returns the square root of this field. Not
   * for use in user code. Use sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = std::sqrt(data);
    return res;
  }

  /**
   * Returns the absolute value of this
   * field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs () const
  {
    VectorizedArray res;
    res.data = std::fabs(data);
    return res;
  }

  /**
   * Returns the component-wise maximum of this
   * field and another one. Not for use in user
   * code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = std::max (data, other.data);
    return res;
  }

  /**
   * Returns the component-wise minimum of this
   * field and another one. Not for use in user
   * code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = std::min (data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2> friend VectorizedArray<Number2>
  std::sqrt (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::abs  (const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::max  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
  template <typename Number2> friend VectorizedArray<Number2>
  std::min  (const VectorizedArray<Number2> &, const VectorizedArray<Number2> &);
};

/**
 * Create a vectorized array that sets all entries in the array to the given
 * scalar.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
make_vectorized_array (const Number &u)
{
  VectorizedArray<Number> result;
  result = u;
  return result;
}

/**
 * Addition of two vectorized arrays with operator +.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator + (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp+=v;
}

/**
 * Subtraction of two vectorized arrays with operator -.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator - (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp-=v;
}

/**
 * Multiplication of two vectorized arrays with operator *.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator * (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp*=v;
}

/**
 * Division of two vectorized arrays with operator /.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator / (const VectorizedArray<Number> &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp = u;
  return tmp/=v;
}

/**
 * Addition of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator + (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp+=v;
}

/**
 * Addition of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator + (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = u;
  return tmp+=v;
}

/**
 * Addition of a vectorized array and a scalar (expanded to a vectorized array
 * with @p n_array_elements equal entries).
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator + (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  return u + v;
}

/**
 * Addition of a vectorized array and a scalar (expanded to a vectorized array
 * with @p n_array_elements equal entries) in case the scalar is a double
 * (needed in order to be able to write simple code with constants that are
 * usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator + (const VectorizedArray<float> &v,
            const double                 &u)
{
  return u + v;
}

/**
 * Subtraction of a vectorized array from a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries).
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator - (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp-=v;
}

/**
 * Subtraction of a vectorized array from a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator - (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp-=v;
}

/**
 * Subtraction of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) from a vectorized array.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator - (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return v-tmp;
}

/**
 * Subtraction of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) from a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator - (const VectorizedArray<float> &v,
            const double                 &u)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return v-tmp;
}

/**
 * Multiplication of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator * (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp*=v;
}

/**
 * Multiplication of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator * (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp*=v;
}

/**
 * Multiplication of a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries).
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator * (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  return u * v;
}

/**
 * Multiplication of a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator * (const VectorizedArray<float> &v,
            const double                 &u)
{
  return u * v;
}

/**
 * Quotient between a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator / (const Number                  &u,
            const VectorizedArray<Number> &v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp/=v;
}

/**
 * Quotient between a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator / (const double                 &u,
            const VectorizedArray<float> &v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp/=v;
}

/**
 * Quotient between a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries).
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator / (const VectorizedArray<Number> &v,
            const Number                  &u)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return v/tmp;
}

/**
 * Quotient between a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relates VectorizedArray
 */
inline
VectorizedArray<float>
operator / (const VectorizedArray<float> &v,
            const double                 &u)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return v/tmp;
}

/**
 * Unary operator + on a vectorized array.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator + (const VectorizedArray<Number> &u)
{
  return u;
}

/**
 * Unary operator - on a vectorized array.
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
VectorizedArray<Number>
operator - (const VectorizedArray<Number> &u)
{
  // to get a negative sign, subtract the input from zero (could also
  // multiply by -1, but this one is slightly simpler)
  return VectorizedArray<Number>()-u;
}



DEAL_II_NAMESPACE_CLOSE


/**
 * Implementation of functions from cmath on VectorizedArray. These functions
 * do not reside in the dealii namespace in order to ensure a similar
 * interface as for the respective functions in cmath. Instead, call them
 * using std::sin.
 */
namespace std
{
  /**
   * Computes the sine of a vectorized data field. The result is returned as
   * vectorized array in the form <tt>{sin(x[0]), sin(x[1]), ...,
   * sin(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  sin (const ::dealii::VectorizedArray<Number> &x)
  {
    // put values in an array and later read in that array with an unaligned
    // read. This should safe some instructions as compared to directly
    // setting the individual elements and also circumvents a compiler
    // optimization bug in gcc-4.6 (see also deal.II developers list from
    // April 2014, topic "matrix_free/step-48 Test").
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      values[i] = std::sin(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * Computes the cosine of a vectorized data field. The result is returned as
   * vectorized array in the form <tt>{cos(x[0]), cos(x[1]), ...,
   * cos(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  cos (const ::dealii::VectorizedArray<Number> &x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      values[i] = std::cos(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * Computes the tangent of a vectorized data field. The result is returned as
   * vectorized array in the form <tt>{tan(x[0]), tan(x[1]), ...,
   * tan(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  tan (const ::dealii::VectorizedArray<Number> &x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      values[i] = std::tan(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * Computes the exponential of a vectorized data field. The result is returned
   * as vectorized array in the form <tt>{exp(x[0]), exp(x[1]), ...,
   * exp(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  exp (const ::dealii::VectorizedArray<Number> &x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      values[i] = std::exp(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * Computes the natural logarithm of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{log(x[0]), log(x[1]), ...,
   * log(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  log (const ::dealii::VectorizedArray<Number> &x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      values[i] = std::log(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }



  /**
   * Computes the square root of a vectorized data field. The result is returned
   * as vectorized array in the form <tt>{sqrt(x[0]), sqrt(x[1]), ...,
   * sqrt(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  sqrt (const ::dealii::VectorizedArray<Number> &x)
  {
    return x.get_sqrt();
  }



  /**
   * Computes the absolute value (modulus) of a vectorized data field. The
   * result is returned as vectorized array in the form <tt>{abs(x[0]),
   * abs(x[1]), ..., abs(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  abs (const ::dealii::VectorizedArray<Number> &x)
  {
    return x.get_abs();
  }



  /**
   * Computes the componentwise maximum of two vectorized data fields. The
   * result is returned as vectorized array in the form <tt>{max(x[0],y[0]),
   * max(x[1],y[1]), ...}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  max (const ::dealii::VectorizedArray<Number> &x,
       const ::dealii::VectorizedArray<Number> &y)
  {
    return x.get_max(y);
  }



  /**
   * Computes the componentwise minimum of two vectorized data fields. The
   * result is returned as vectorized array in the form <tt>{min(x[0],y[0]),
   * min(x[1],y[1]), ...}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  min (const ::dealii::VectorizedArray<Number> &x,
       const ::dealii::VectorizedArray<Number> &y)
  {
    return x.get_min(y);
  }

}

#endif
