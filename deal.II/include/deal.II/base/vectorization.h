// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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
#include <deal.II/base/std_cxx1x/type_traits.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/parallel.h>

#include <cmath>
#include <cstring>

// Note:
// The flag DEAL_II_COMPILER_VECTORIZATION_LEVEL is essentially constructed
// according to the following scheme
// #ifdef __AVX__
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 2
// #elif defined (__SSE2__)
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 1
// #else
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 0
// #endif
// In addition to checking the flags __AVX__ and __SSE2__, a configure test
// ensures that these feature are not only present but also working properly.

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL == 2 // AVX
#include <immintrin.h>
#include <mm_malloc.h>
#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL == 1 // SSE2
#include <emmintrin.h>
#include <mm_malloc.h>
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


// for safety, also check that __AVX__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL == 2  && defined(__AVX__)

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



/**
 * This namespace defines the copy and set functions used in
 * AlignedVector. These functions operate in parallel when there are enough
 * elements in the vector.
 */
namespace internal
{
  /**
   * Move and class that actually issues the copy commands in
   * AlignedVector. This class is based on the specialized for loop base class
   * ParallelForLoop in parallel.h whose purpose is the following: When
   * calling a parallel for loop on AlignedVector with apply_to_subranges, it
   * generates different code for every different argument we might choose (as
   * it is templated). This gives a lot of code (e.g. it triples the memory
   * required for compiling the file matrix_free.cc and the final object size
   * is several times larger) which is completely useless. Therefore, this
   * class channels all copy commands through one call to apply_to_subrange
   * for all possible types, which makes the copy operation much cleaner
   * (thanks to a virtual function, whose cost is negligible in this context).
   *
   * @relates AlignedVector
   */
  template <typename T>
  class AlignedVectorMove : private parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size = 160000/sizeof(T)+1;
  public:
    /**
     * Constructor. Issues a parallel call if
     * there are sufficiently many elements,
     * otherwise work in serial. Copies the data
     * from source to destination and then calls
     * destructor on the source. If the optional
     * argument is set to true, the source is left
     * untouched instead.
     */
    AlignedVectorMove (T *source_begin,
                       T *source_end,
                       T *destination,
                       bool copy_only = false)
      :
      source_ (source_begin),
      destination_ (destination),
      copy_only_ (copy_only)
    {
      Assert (source_end >= source_begin, ExcInternalError());
      const std::size_t size = source_end - source_begin;
      if (size < minimum_parallel_grain_size)
        apply_to_subrange (0, size);
      else
        apply_parallel (0, size, minimum_parallel_grain_size);
    }

    /**
     * This method moves elements from the source
     * to the destination given in the constructor
     * on a subrange given by two integers.
     */
    virtual void apply_to_subrange (const std::size_t begin,
                                    const std::size_t end) const
    {
      // for classes trivial assignment can use
      // memcpy
      if (std_cxx1x::is_trivial<T>::value == true)
        std::memcpy (destination_+begin, source_+begin, (end-begin)*sizeof(T));
      else if (copy_only_ == false)
        for (std::size_t i=begin; i<end; ++i)
          {
            // initialize memory, copy, and destruct
            new (&destination_[i]) T;
            destination_[i] = source_[i];
            source_[i].~T();
          }
      else
        for (std::size_t i=begin; i<end; ++i)
          {
            new (&destination_[i]) T;
            destination_[i] = source_[i];
          }
    }

  private:
    T *source_;
    T *destination_;
    const bool copy_only_;
  };

  /**
   * Class that issues the set commands for AlignedVector.
   *
   * @relates AlignedVector
   */
  template <typename T>
  class AlignedVectorSet : private parallel::ParallelForInteger
  {
    static const std::size_t minimum_parallel_grain_size = 160000/sizeof(T)+1;
  public:
    /**
     * Constructor. Issues a parallel call if
     * there are sufficiently many elements,
     * otherwise work in serial.
     */
    AlignedVectorSet (const std::size_t size,
                      const T &element,
                      T *destination)
      :
      element_ (element),
      destination_ (destination),
      trivial_element (false)
    {
      if (size == 0)
        return;

      if (std_cxx1x::is_trivial<T>::value == true)
        {
          const unsigned char zero [sizeof(T)] = {};
          if (std::memcmp(zero, &element, sizeof(T)) == 0)
            trivial_element = true;
        }
      if (size < minimum_parallel_grain_size)
        apply_to_subrange (0, size);
      else
        apply_parallel (0, size, minimum_parallel_grain_size);
    }

  private:

    /**
     * This sets elements on a subrange given by
     * two integers.
     */
    virtual void apply_to_subrange (const std::size_t begin,
                                    const std::size_t end) const
    {
      // for classes with trivial assignment of zero
      // can use memset
      if (std_cxx1x::is_trivial<T>::value == true && trivial_element)
        std::memset (destination_+begin, 0, (end-begin)*sizeof(T));
      else
        for (std::size_t i=begin; i<end; ++i)
          {
            // initialize memory and set
            new (&destination_[i]) T;
            destination_[i] = element_;
          }
    }

    const T &element_;
    mutable T *destination_;
    bool trivial_element;
  };
} // end of namespace internal


/**
 * This is a replacement class for std::vector to be used in combination with
 * VectorizedArray and derived data types. It allocates memory aligned to
 * addresses of a vectorized data type (for SSE, this is necessary in order to
 * avoid segfaults, and for AVX it considerably increases performance). This
 * could also be achieved by proving std::vector with a user-defined
 * allocator. On the other hand, writing an own small vector class lets us
 * insert assertions more easily, and cut some unnecessary functionality. Note
 * that this vector is a bit more memory-consuming than std::vector because of
 * alignment, so it is recommended to only use this vector on long vectors.
 *
 * @p author Katharina Kormann, Martin Kronbichler, 2011
 */
template < class T >
class AlignedVector
{
public:
  /**
   * Declare standard types used in all
   * containers. These types parallel those
   * in the <tt>C++</tt> standard libraries
   * <tt>vector<...></tt> class.
   */
  typedef T                   value_type;
  typedef value_type         *pointer;
  typedef const value_type   *const_pointer;
  typedef value_type         *iterator;
  typedef const value_type   *const_iterator;
  typedef value_type         &reference;
  typedef const value_type   &const_reference;
  typedef std::size_t         size_type;

  /**
   * Empty constructor. Sets the vector size to
   * zero.
   */
  AlignedVector ()
    :
    _data (0),
    _end_data (0),
    _end_allocated (0)
  {};

  /**
   * Sets the vector size to the given size and
   * initializes all elements with T().
   */
  AlignedVector (const size_type size)
    :
    _data (0),
    _end_data (0),
    _end_allocated (0)
  {
    if (size > 0)
      resize (size);
  }

  /**
   * Destructor.
   */
  ~AlignedVector ()
  {
    clear();
  }

  /**
   * Copy constructor.
   */
  AlignedVector (const AlignedVector<T> &vec)
    :
    _data (0),
    _end_data (0),
    _end_allocated (0)
  {
    // do not invalidate old data
    resize_fast (vec._end_data - vec._data);
    internal::AlignedVectorMove<T> (vec._data, vec._end_data, _data, true);
  }

  /**
   * Assignment to the input vector @p vec.
   */
  AlignedVector &
  operator = (const AlignedVector<T> &vec)
  {
    clear();
    resize_fast (vec._end_data - vec._data);
    internal::AlignedVectorMove<T> (vec._data, vec._end_data, _data, true);
    return *this;
  }

  /**
   * Change the size of the vector. It keeps old
   * elements previously available but does not
   * initialize the newly allocated memory,
   * leaving it in an undefined state.
   */
  void resize_fast (const size_type size)
  {
    reserve (size);
    _end_data = _data + size;
  }

  /**
   * Change the size of the vector. It keeps old
   * elements previously available, and
   * initializes each element with the specified
   * data. If the new vector size is shorter
   * than the old one, the memory is not
   * released unless the new size is zero.
   */
  void resize (const size_type size_in,
               const T        &init = T())
  {
    const size_type old_size = size();
    if (std_cxx1x::is_trivial<T>::value == false && size_in < old_size)
      {
        // call destructor on fields that are released
        while (_end_data != _data+size_in)
          (--_end_data)->~T();
      }

    resize_fast (size_in);
    // now _size is set correctly, need to set the
    // values
    if (size_in > old_size)
      internal::AlignedVectorSet<T> (size_in-old_size, init,
                                     _data+old_size);
  }

  /**
   * Reserve memory space for @p size
   * elements. If the argument @p size is set to
   * zero, all previously allocated memory is
   * released.
   *
   * In order to avoid too frequent reallocation
   * (which involves copy of the data), this
   * function doubles the amount of memory
   * occupied when the given size is larger than
   * the previously allocated size.
   */
  void reserve (const size_type size_alloc)
  {
    const size_type old_size = _end_data - _data;
    const size_type allocated_size = _end_allocated - _data;
    if (size_alloc > allocated_size)
      {
        // if we continuously increase the size of the
        // vector, we might be reallocating a lot of
        // times. therefore, try to increase the size
        // more aggressively
        size_type new_size = size_alloc;
        if (size_alloc < (2 * allocated_size))
          new_size = 2 * allocated_size;

        const size_type size_actual_allocate = new_size * sizeof(T);

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL > 0

        // allocate and align along boundaries of the
        // size of VectorizedArray<double>, which is
        // 16 bytes for SSE and 32 bytes for AVX
        T *new_data = static_cast<T *>(_mm_malloc (size_actual_allocate,
                                                   sizeof(VectorizedArray<double>)));
#else
        T *new_data = static_cast<T *>(malloc (size_actual_allocate));
#endif
        if (new_data == 0)
          throw std::bad_alloc();

        // copy data in case there was some content
        // before and release the old memory with the
        // function corresponding to the one used for
        // allocating
        std::swap (_data, new_data);
        _end_data = _data + old_size;
        _end_allocated = _data + new_size;
        if (_end_data != _data)
          {
            internal::AlignedVectorMove<T>(new_data, new_data + old_size,
                                           _data);
#if DEAL_II_COMPILER_VECTORIZATION_LEVEL > 0
            _mm_free(new_data);
#else
            free(new_data);
#endif
          }
      }
    else if (size_alloc == 0)
      clear();
  }

  /**
   * Releases all previously allocated memory
   * and leaves the vector in a state equivalent
   * to the state after the default constructor
   * has been called.
   */
  void clear ()
  {
    if (_data != 0)
      {
        if (std_cxx1x::is_trivial<T>::value == false)
          while (_end_data != _data)
            (--_end_data)->~T();

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL > 0
        _mm_free(_data);
#else
        free(_data);
#endif
      }
    _data = 0;
    _end_data = 0;
    _end_allocated = 0;
  };

  /**
   * Inserts an element at the end of the
   * vector, increasing the vector size by
   * one. Note that the allocated size will
   * double whenever the previous space is not
   * enough to hold the new element.
   */
  void push_back (const T in_data)
  {
    Assert (_end_data <= _end_allocated, ExcInternalError());
    if (_end_data == _end_allocated)
      reserve (std::max(2*capacity(),static_cast<size_type>(16)));
    if (std_cxx1x::is_trivial<T>::value == false)
      new (_end_data) T;
    *_end_data++ = in_data;
  }

  /**
   * Returns the last element of the vector
   * (read and write access).
   */
  reference back ()
  {
    AssertIndexRange (0, size());
    T *field = _end_data - 1;
    return *field;
  }

  /**
   * Returns the last element of the vector
   * (read-only access).
   */
  const_reference back () const
  {
    AssertIndexRange (0, size());
    const T *field = _end_data - 1;
    return *field;
  }

  /**
   * Inserts several elements at the end of the
   * vector given by a range of elements.
   */
  template <typename ForwardIterator>
  void insert_back (ForwardIterator begin,
                    ForwardIterator end)
  {
    const unsigned int old_size = size();
    reserve (old_size + (end-begin));
    for ( ; begin != end; ++begin, ++_end_data)
      {
        if (std_cxx1x::is_trivial<T>::value == false)
          new (_end_data) T;
        *_end_data = *begin;
      }
  }

  /**
   * Swaps the given vector with the calling
   * vector.
   */
  void swap (AlignedVector<T> &vec)
  {
    std::swap (_data, vec._data);
    std::swap (_end_data, vec._end_data);
    std::swap (_end_allocated, vec._end_allocated);
  }

  /**
   * Returns the size of the vector.
   */
  size_type size () const
  {
    return _end_data - _data;
  }

  /**
   * Returns the capacity of the vector, i.e.,
   * the size this vector can hold without
   * reallocation. Note that capacity() >=
   * size().
   */
  size_type capacity () const
  {
    return _end_allocated - _data;
  }

  /**
   * Read-write access to entry @p index in the
   * vector.
   */
  reference
  operator [] (const size_type index)
  {
    AssertIndexRange (index, size());
    return _data[index];
  };

  /**
   * Read-only access to entry @p index in the
   * vector.
   */
  const_reference operator [] (const size_type index) const
  {
    AssertIndexRange (index, size());
    return _data[index];
  };

  /**
   * Returns a read and write pointer to the
   * beginning of the data array.
   */
  iterator begin ()
  {
    return _data;
  }

  /**
   * Returns a read and write pointer to the
   * end of the data array.
   */
  iterator end ()
  {
    return _end_data;
  }

  /**
   * Returns a read-only pointer to the
   * beginning of the data array.
   */
  const_iterator begin () const
  {
    return _data;
  }

  /**
   * Returns a read-only pointer to the
   * end of the data array.
   */
  const_iterator end () const
  {
    return _end_data;
  }

  /**
   * Returns the memory consumption of the
   * allocated memory in this class. If the
   * underlying type @p T allocates memory by
   * itself, this memory is not counted.
   */
  size_type memory_consumption () const
  {
    size_type memory = sizeof(this);
    memory += sizeof(T) * capacity();
    return memory;
  }

private:

  /**
   * Pointer to actual class data.
   */
  T *_data;

  /**
   * Pointer to the end of valid data fields.
   */
  T *_end_data;

  /**
   * Pointer to the end of the allocated memory.
   */
  T *_end_allocated;
};


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
   * Computes the sine of a vectorized
   * data field. The result is return as
   * vectorized array in the form
   * <tt>{sin(x[0]), sin(x[1]), ...,
   * sin(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  sin (const ::dealii::VectorizedArray<Number> &x)
  {
    ::dealii::VectorizedArray<Number> sin_val;
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      sin_val[i] = std::sin(x[i]);
    return sin_val;
  }



  /**
   * Computes the tangent of a vectorized
   * data field. The result is return as
   * vectorized array in the form
   * <tt>{tan(x[0]), tan(x[1]), ...,
   * tan(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  tan (const ::dealii::VectorizedArray<Number> &x)
  {
    ::dealii::VectorizedArray<Number> tan_val;
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      tan_val[i] = std::tan(x[i]);
    return tan_val;
  }


  /**
   * Computes the cosine of a vectorized
   * data field. The result is return as
   * vectorized array in the form
   * <tt>{cos(x[0]), cos(x[1]), ...,
   * cos(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  cos (const ::dealii::VectorizedArray<Number> &x)
  {
    ::dealii::VectorizedArray<Number> cos_val;
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      cos_val[i] = std::cos(x[i]);
    return cos_val;
  }


  /**
   * Computes the exponential of a vectorized
   * data field. The result is return as
   * vectorized array in the form
   * <tt>{exp(x[0]), exp(x[1]), ...,
   * exp(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  exp (const ::dealii::VectorizedArray<Number> &x)
  {
    ::dealii::VectorizedArray<Number> exp_val;
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      exp_val[i] = std::exp(x[i]);
    return exp_val;
  }


  /**
   * Computes the natural logarithm of a
   * vectorized data field. The result is return
   * as vectorized array in the form
   * <tt>{log(x[0]), log(x[1]), ...,
   * log(x[n_array_elements-1])}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  log (const ::dealii::VectorizedArray<Number> &x)
  {
    ::dealii::VectorizedArray<Number> log_val;
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      log_val[i] = std::log(x[i]);
    return log_val;
  }



  /**
   * Computes the square root of a vectorized
   * data field. The result is return as
   * vectorized array in the form
   * <tt>{sqrt(x[0]), sqrt(x[1]), ...,
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
   * Computes the absolute value (modulus) of a
   * vectorized data field. The result is return
   * as vectorized array in the form
   * <tt>{abs(x[0]), abs(x[1]), ...,
   * abs(x[n_array_elements-1])}</tt>.
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
   * Computes the componentwise maximum of two
   * vectorized data fields. The result is
   * return as vectorized array in the form
   * <tt>{max(x[0],y[0]), max(x[1],y[1]),
   * ...}</tt>.
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
   * Computes the componentwise minimum of two
   * vectorized data fields. The result is
   * return as vectorized array in the form
   * <tt>{min(x[0],y[0]), min(x[1],y[1]),
   * ...}</tt>.
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
