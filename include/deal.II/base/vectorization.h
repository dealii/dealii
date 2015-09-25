// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
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


#ifndef dealii__vectorization_h
#define dealii__vectorization_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

#include <cmath>

// Note:
// The flag DEAL_II_COMPILER_VECTORIZATION_LEVEL is essentially constructed
// according to the following scheme
// #ifdef __AVX512F__
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 3
// #elif defined (__AVX__)
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 2
// #elif defined (__SSE2__)
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 1
// #else
// #define DEAL_II_COMPILER_VECTORIZATION_LEVEL 0
// #endif
// In addition to checking the flags __AVX__ and __SSE2__, a CMake test,
// 'check_01_cpu_features.cmake', ensures that these feature are not only
// present in the compilation unit but also working properly.

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 // AVX, AVX-512
#include <immintrin.h>
#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL == 1 // SSE2
#include <emmintrin.h>
#endif


// forward declarations
DEAL_II_NAMESPACE_OPEN
template <typename Number> class VectorizedArray;
template <typename T> struct EnableIfScalar;
DEAL_II_NAMESPACE_CLOSE


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


// Enable the EnableIfScalar type trait for VectorizedArray<Number> such
// that it can be used as a Number type in Tensor<rank,dim,Number>, etc.

template<typename Number>
struct EnableIfScalar<VectorizedArray<Number> >
{
  typedef VectorizedArray<typename EnableIfScalar<Number>::type> type;
};


/**
 * This generic class defines a unified interface to a vectorized data type.
 * For general template arguments, this class simply corresponds to the
 * template argument. For example, VectorizedArray<long double> is nothing
 * else but a wrapper around <tt>long double</tt> with exactly one data field
 * of type <tt>long double</tt> and overloaded arithmetic operations. This
 * means that <tt>VectorizedArray<ComplicatedType></tt> has a similar layout
 * as ComplicatedType, provided that ComplicatedType defines basic arithmetic
 * operations. For floats and doubles, an array of numbers are packed
 * together, though. The number of elements packed together depend on the
 * computer system and compiler flags that are used for compilation of
 * deal.II. The fundamental idea of these packed data types is to use one
 * single CPU instruction to perform arithmetic operations on the whole array
 * using the processor's vector units. Most computer systems by 2010 standards
 * will use an array of two doubles and four floats, respectively (this
 * corresponds to the SSE/SSE2 data sets) when compiling deal.II on 64-bit
 * operating systems. On Intel Sandy Bridge processors and newer or AMD
 * Bulldozer processors and newer, four doubles and eight floats are used when
 * deal.II is configured e.g. using gcc with --with-cpu=native or --with-
 * cpu=corei7-avx. On compilations with AVX-512 support, eight doubles and
 * sixteen floats are used.
 *
 * This behavior of this class is made similar to the basic data types double
 * and float. The definition of a vectorized array does not initialize the
 * data field but rather leaves it undefined, as is the case for double and
 * float. However, when calling something like VectorizedArray<double> a =
 * VectorizedArray<double>(), it sets all numbers in this field to zero. In
 * other words, this class is a plain old data (POD) type which has an
 * equivalent C representation and can e.g. be safely copied with std::memcpy.
 * This POD layout is also necessary for ensuring correct alignment of data
 * with address boundaries when collected in a vector (i.e., when the first
 * element in a vector is properly aligned, all subsequent elements will be
 * correctly aligned, too).
 *
 * Note that for proper functioning of this class, certain data alignment
 * rules must be respected. This is because the computer expects the starting
 * address of a VectorizedArray<double> field at specific addresses in memory
 * (usually, the address of the vectorized array should be a multiple of the
 * length of the array in bytes). Otherwise, a segmentation fault or a severe
 * loss of performance might occur. When creating a single data field on the
 * stack like <tt>VectorizedArray<double> a = VectorizedArray<double>()</tt>,
 * the compiler will take care of data alignment automatically. However, when
 * allocating a long vector of VectorizedArray<double> data, one needs to
 * respect these rules. Use the class AlignedVector or data containers based
 * on AlignedVector (such as Table) for this purpose. It is a class very
 * similar to std::vector otherwise but always makes sure that data is
 * correctly aligned.
 *
 * @author Katharina Kormann, Martin Kronbichler, 2010, 2011
 */
template <typename Number>
class VectorizedArray
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 1;

  // POD means that there should be no user-defined constructors, destructors
  // and copy functions (the standard is somewhat relaxed in C++2011, though).

  /**
   * This function assigns a scalar to this class.
   */
  VectorizedArray &
  operator = (const Number scalar)
  {
    data = scalar;
    return *this;
  }

  /**
   * Access operator (only valid with component 0)
   */
  Number &
  operator [] (const unsigned int comp)
  {
    (void)comp;
    AssertIndexRange (comp, 1);
    return data;
  }

  /**
   * Constant access operator (only valid with component 0)
   */
  const Number &
  operator [] (const unsigned int comp) const
  {
    (void)comp;
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
   * Actual data field. Since this class represents a POD data type, it is
   * declared public.
   */
  Number data;

private:
  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = std::sqrt(data);
    return res;
  }

  /**
   * Returns the absolute value of this field. Not for use in user code. Use
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
   * Returns the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max (const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = std::max (data, other.data);
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
 * This method loads VectorizedArray::n_array_elements data streams from the
 * given array @p in. The offsets to the input array are given by the array @p
 * offsets. From each stream, n_entries are read. The data is then transposed
 * and stored it into an array of VectorizedArray type. The output array @p
 * out is expected to be an array of size @p n_entries. This method operates
 * on plain arrays, so no checks for valid data access are made. It is the
 * user's responsibility to ensure that the given arrays are valid according
 * to the access layout below.
 *
 * This operation corresponds to a transformation of an array-of-struct
 * (input) into a struct-of-array (output) according to the following formula:
 *
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 *   for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
 *     out[i][v] = in[offsets[v]+i];
 * @endcode
 *
 * A more optimized version of this code will be used for supported types.
 *
 * This is the inverse operation to vectorized_transpose_and_store().
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const Number            *in,
                              const unsigned int      *offsets,
                              VectorizedArray<Number> *out)
{
  for (unsigned int i=0; i<n_entries; ++i)
    for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
      out[i][v] = in[offsets[v]+i];
}



/**
 * This method stores the vectorized arrays in transposed form into the given
 * output array @p out with the given offsets @p offsets. This operation
 * corresponds to a transformation of a struct-of-array (input) into an array-
 * of-struct (output). This method operates on plain array, so no checks for
 * valid data access are made. It is the user's responsibility to ensure that
 * the given arrays are valid according to the access layout below.
 *
 * This method assumes that the specified offsets do not overlap. Otherwise,
 * the behavior is undefined in the vectorized case. It is the user's
 * responsibility to make sure that the access does not overlap and avoid
 * undefined behavior.
 *
 * The argument @p add_into selects where the entries should only be written
 * into the output arrays or the result should be added into the existing
 * entries in the output. For <code>add_into == false</code>, the following
 * code is assumed:
 *
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 *   for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
 *     out[offsets[v]+i] = in[i][v];
 * @endcode
 *
 * For <code>add_into == true</code>, the code implements the following
 * action:
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 *   for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
 *     out[offsets[v]+i] += in[i][v];
 * @endcode
 *
 * A more optimized version of this code will be used for supported types.
 *
 * This is the inverse operation to vectorized_load_and_transpose().
 *
 * @relates VectorizedArray
 */
template <typename Number>
inline
void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<Number> *in,
                               const unsigned int            *offsets,
                               Number                        *out)
{
  if (add_into)
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        out[offsets[v]+i] = in[i][v];
}



// for safety, also check that __AVX512F__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 3  && defined(__AVX512F__)

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
    data = _mm512_set1_pd(x);
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
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
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
    // value to +. Since there is no andnot for AVX512, we interpret the data
    // as 64 bit integers and do the andnot on those types (note that andnot
    // is a bitwise operation so the data type does not matter)
    __m512d mask = _mm512_set1_pd (-0.);
    VectorizedArray res;
    res.data = (__m512d)_mm512_andnot_epi64 ((__m512i)mask, (__m512i)data);
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
    data = _mm512_set1_ps(x);
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
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
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
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +. Since there is no andnot for AVX512, we interpret the data
    // as 32 bit integers and do the andnot on those types (note that andnot
    // is a bitwise operation so the data type does not matter)
    __m512 mask = _mm512_set1_ps (-0.f);
    VectorizedArray res;
    res.data = (__m512)_mm512_andnot_epi32 ((__m512i)mask, (__m512i)data);
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



#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2  && defined(__AVX__)

/**
 * Specialization of VectorizedArray class for double and AVX.
 */
template <>
class VectorizedArray<double>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 4;

  /**
   * This function can be used to set all data fields to a given scalar.
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
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
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
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m256d data;

private:
  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_pd(data);
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
    __m256d mask = _mm256_set1_pd (-0.);
    VectorizedArray res;
    res.data = _mm256_andnot_pd(mask, data);
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
    res.data = _mm256_max_pd (data, other.data);
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
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const double            *in,
                              const unsigned int      *offsets,
                              VectorizedArray<double> *out)
{
  const unsigned int n_chunks = n_entries/4, remainder = n_entries%4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m256d u0 = _mm256_loadu_pd(in+4*i+offsets[0]);
      __m256d u1 = _mm256_loadu_pd(in+4*i+offsets[1]);
      __m256d u2 = _mm256_loadu_pd(in+4*i+offsets[2]);
      __m256d u3 = _mm256_loadu_pd(in+4*i+offsets[3]);
      __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
      __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
      __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
      __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
      out[4*i+0].data = _mm256_unpacklo_pd (t0, t1);
      out[4*i+1].data = _mm256_unpackhi_pd (t0, t1);
      out[4*i+2].data = _mm256_unpacklo_pd (t2, t3);
      out[4*i+3].data = _mm256_unpackhi_pd (t2, t3);
    }
  if (remainder > 0 && n_chunks > 0)
    {
      // simple re-load all data in the last slot
      const unsigned int final_pos = n_chunks*4-4+remainder;
      Assert(final_pos+4 == n_entries, ExcInternalError());
      __m256d u0 = _mm256_loadu_pd(in+final_pos+offsets[0]);
      __m256d u1 = _mm256_loadu_pd(in+final_pos+offsets[1]);
      __m256d u2 = _mm256_loadu_pd(in+final_pos+offsets[2]);
      __m256d u3 = _mm256_loadu_pd(in+final_pos+offsets[3]);
      __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
      __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
      __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
      __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
      out[final_pos+0].data = _mm256_unpacklo_pd (t0, t1);
      out[final_pos+1].data = _mm256_unpackhi_pd (t0, t1);
      out[final_pos+2].data = _mm256_unpacklo_pd (t2, t3);
      out[final_pos+3].data = _mm256_unpackhi_pd (t2, t3);
    }
  else if (remainder > 0)
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[i][v] = in[offsets[v]+i];
}



/**
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double> *in,
                               const unsigned int            *offsets,
                               double                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m256d u0 = in[4*i+0].data;
      __m256d u1 = in[4*i+1].data;
      __m256d u2 = in[4*i+2].data;
      __m256d u3 = in[4*i+3].data;
      __m256d t0 = _mm256_permute2f128_pd (u0, u2, 0x20);
      __m256d t1 = _mm256_permute2f128_pd (u1, u3, 0x20);
      __m256d t2 = _mm256_permute2f128_pd (u0, u2, 0x31);
      __m256d t3 = _mm256_permute2f128_pd (u1, u3, 0x31);
      __m256d res0 = _mm256_unpacklo_pd (t0, t1);
      __m256d res1 = _mm256_unpackhi_pd (t0, t1);
      __m256d res2 = _mm256_unpacklo_pd (t2, t3);
      __m256d res3 = _mm256_unpackhi_pd (t2, t3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out+4*i+offsets[0]), res0);
          _mm256_storeu_pd(out+4*i+offsets[0], res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out+4*i+offsets[1]), res1);
          _mm256_storeu_pd(out+4*i+offsets[1], res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out+4*i+offsets[2]), res2);
          _mm256_storeu_pd(out+4*i+offsets[2], res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out+4*i+offsets[3]), res3);
          _mm256_storeu_pd(out+4*i+offsets[3], res3);
        }
      else
        {
          _mm256_storeu_pd(out+4*i+offsets[0], res0);
          _mm256_storeu_pd(out+4*i+offsets[1], res1);
          _mm256_storeu_pd(out+4*i+offsets[2], res2);
          _mm256_storeu_pd(out+4*i+offsets[3], res3);
        }
    }
  const unsigned int shift = n_chunks * 4;
  if (add_into)
    for (unsigned int i=shift; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=shift; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[offsets[v]+i] = in[i][v];
}



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
  static const unsigned int n_array_elements = 8;

  /**
   * This function can be used to set all data fields to a given scalar.
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
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
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
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m256 data;

private:

  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_ps(data);
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
    __m256 mask = _mm256_set1_ps (-0.f);
    VectorizedArray res;
    res.data = _mm256_andnot_ps(mask, data);
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
    res.data = _mm256_max_ps (data, other.data);
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



/**
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_load_and_transpose(const unsigned int      n_entries,
                              const float            *in,
                              const unsigned int     *offsets,
                              VectorizedArray<float> *out)
{
  const unsigned int n_chunks = n_entries/4, remainder = n_entries%4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m128 u0 = _mm_loadu_ps(in+4*i+offsets[0]);
      __m128 u1 = _mm_loadu_ps(in+4*i+offsets[1]);
      __m128 u2 = _mm_loadu_ps(in+4*i+offsets[2]);
      __m128 u3 = _mm_loadu_ps(in+4*i+offsets[3]);
      __m128 u4 = _mm_loadu_ps(in+4*i+offsets[4]);
      __m128 u5 = _mm_loadu_ps(in+4*i+offsets[5]);
      __m128 u6 = _mm_loadu_ps(in+4*i+offsets[6]);
      __m128 u7 = _mm_loadu_ps(in+4*i+offsets[7]);
      // To avoid warnings about uninitialized variables, need to initialize
      // one variable with zero before using it.
      __m256 t0, t1, t2, t3 = _mm256_set1_ps(0.F);
      t0 = _mm256_insertf128_ps (t3, u0, 0);
      t0 = _mm256_insertf128_ps (t0, u4, 1);
      t1 = _mm256_insertf128_ps (t3, u1, 0);
      t1 = _mm256_insertf128_ps (t1, u5, 1);
      t2 = _mm256_insertf128_ps (t3, u2, 0);
      t2 = _mm256_insertf128_ps (t2, u6, 1);
      t3 = _mm256_insertf128_ps (t3, u3, 0);
      t3 = _mm256_insertf128_ps (t3, u7, 1);
      __m256 v0 = _mm256_shuffle_ps (t0, t1, 0x44);
      __m256 v1 = _mm256_shuffle_ps (t0, t1, 0xee);
      __m256 v2 = _mm256_shuffle_ps (t2, t3, 0x44);
      __m256 v3 = _mm256_shuffle_ps (t2, t3, 0xee);
      out[4*i+0].data = _mm256_shuffle_ps (v0, v2, 0x88);
      out[4*i+1].data = _mm256_shuffle_ps (v0, v2, 0xdd);
      out[4*i+2].data = _mm256_shuffle_ps (v1, v3, 0x88);
      out[4*i+3].data = _mm256_shuffle_ps (v1, v3, 0xdd);
    }
  if (remainder > 0 && n_chunks > 0)
    {
      // simple re-load all data in the last slot
      const unsigned int final_pos = n_chunks*4-4+remainder;
      Assert(final_pos+4 == n_entries, ExcInternalError());
      __m128 u0 = _mm_loadu_ps(in+final_pos+offsets[0]);
      __m128 u1 = _mm_loadu_ps(in+final_pos+offsets[1]);
      __m128 u2 = _mm_loadu_ps(in+final_pos+offsets[2]);
      __m128 u3 = _mm_loadu_ps(in+final_pos+offsets[3]);
      __m128 u4 = _mm_loadu_ps(in+final_pos+offsets[4]);
      __m128 u5 = _mm_loadu_ps(in+final_pos+offsets[5]);
      __m128 u6 = _mm_loadu_ps(in+final_pos+offsets[6]);
      __m128 u7 = _mm_loadu_ps(in+final_pos+offsets[7]);
      __m256 t0, t1, t2, t3 = _mm256_set1_ps(0.F);
      t0 = _mm256_insertf128_ps (t3, u0, 0);
      t0 = _mm256_insertf128_ps (t0, u4, 1);
      t1 = _mm256_insertf128_ps (t3, u1, 0);
      t1 = _mm256_insertf128_ps (t1, u5, 1);
      t2 = _mm256_insertf128_ps (t3, u2, 0);
      t2 = _mm256_insertf128_ps (t2, u6, 1);
      t3 = _mm256_insertf128_ps (t3, u3, 0);
      t3 = _mm256_insertf128_ps (t3, u7, 1);
      __m256 v0 = _mm256_shuffle_ps (t0, t1, 0x44);
      __m256 v1 = _mm256_shuffle_ps (t0, t1, 0xee);
      __m256 v2 = _mm256_shuffle_ps (t2, t3, 0x44);
      __m256 v3 = _mm256_shuffle_ps (t2, t3, 0xee);
      out[final_pos+0].data = _mm256_shuffle_ps (v0, v2, 0x88);
      out[final_pos+1].data = _mm256_shuffle_ps (v0, v2, 0xdd);
      out[final_pos+2].data = _mm256_shuffle_ps (v1, v3, 0x88);
      out[final_pos+3].data = _mm256_shuffle_ps (v1, v3, 0xdd);
    }
  else if (remainder > 0)
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<8; ++v)
        out[i][v] = in[offsets[v]+i];
}



/**
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float> *in,
                               const unsigned int           *offsets,
                               float                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m256 u0 = in[4*i+0].data;
      __m256 u1 = in[4*i+1].data;
      __m256 u2 = in[4*i+2].data;
      __m256 u3 = in[4*i+3].data;
      __m256 t0 = _mm256_shuffle_ps (u0, u1, 0x44);
      __m256 t1 = _mm256_shuffle_ps (u0, u1, 0xee);
      __m256 t2 = _mm256_shuffle_ps (u2, u3, 0x44);
      __m256 t3 = _mm256_shuffle_ps (u2, u3, 0xee);
      u0 = _mm256_shuffle_ps (t0, t2, 0x88);
      u1 = _mm256_shuffle_ps (t0, t2, 0xdd);
      u2 = _mm256_shuffle_ps (t1, t3, 0x88);
      u3 = _mm256_shuffle_ps (t1, t3, 0xdd);
      __m128 res0 = _mm256_extractf128_ps (u0, 0);
      __m128 res4 = _mm256_extractf128_ps (u0, 1);
      __m128 res1 = _mm256_extractf128_ps (u1, 0);
      __m128 res5 = _mm256_extractf128_ps (u1, 1);
      __m128 res2 = _mm256_extractf128_ps (u2, 0);
      __m128 res6 = _mm256_extractf128_ps (u2, 1);
      __m128 res3 = _mm256_extractf128_ps (u3, 0);
      __m128 res7 = _mm256_extractf128_ps (u3, 1);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[0]), res0);
          _mm_storeu_ps(out+4*i+offsets[0], res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[1]), res1);
          _mm_storeu_ps(out+4*i+offsets[1], res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[2]), res2);
          _mm_storeu_ps(out+4*i+offsets[2], res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[3]), res3);
          _mm_storeu_ps(out+4*i+offsets[3], res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[4]), res4);
          _mm_storeu_ps(out+4*i+offsets[4], res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[5]), res5);
          _mm_storeu_ps(out+4*i+offsets[5], res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[6]), res6);
          _mm_storeu_ps(out+4*i+offsets[6], res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[7]), res7);
          _mm_storeu_ps(out+4*i+offsets[7], res7);
        }
      else
        {
          _mm_storeu_ps(out+4*i+offsets[0], res0);
          _mm_storeu_ps(out+4*i+offsets[1], res1);
          _mm_storeu_ps(out+4*i+offsets[2], res2);
          _mm_storeu_ps(out+4*i+offsets[3], res3);
          _mm_storeu_ps(out+4*i+offsets[4], res4);
          _mm_storeu_ps(out+4*i+offsets[5], res5);
          _mm_storeu_ps(out+4*i+offsets[6], res6);
          _mm_storeu_ps(out+4*i+offsets[7], res7);
        }
    }
  const unsigned int shift = n_chunks * 4;
  if (add_into)
    for (unsigned int i=shift; i<n_entries; ++i)
      for (unsigned int v=0; v<8; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=shift; i<n_entries; ++i)
      for (unsigned int v=0; v<8; ++v)
        out[offsets[v]+i] = in[i][v];
}



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
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 2;

  /**
   * This function can be used to set all data fields to a given scalar.
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
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m128d data;

private:
  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_pd(data);
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
    __m128d mask = _mm_set1_pd (-0.);
    VectorizedArray res;
    res.data = _mm_andnot_pd(mask, data);
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
    res.data = _mm_max_pd (data, other.data);
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
 * Specialization for double and SSE2.
 */
template <>
inline
void vectorized_load_and_transpose(const unsigned int      n_entries,
                                   const double            *in,
                                   const unsigned int      *offsets,
                                   VectorizedArray<double> *out)
{
  const unsigned int n_chunks = n_entries/2, remainder = n_entries%2;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m128d u0 = _mm_loadu_pd(in+2*i+offsets[0]);
      __m128d u1 = _mm_loadu_pd(in+2*i+offsets[1]);
      out[2*i+0].data = _mm_unpacklo_pd (u0, u1);
      out[2*i+1].data = _mm_unpackhi_pd (u0, u1);
    }
  if (remainder > 0)
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<2; ++v)
        out[i][v] = in[offsets[v]+i];
}



/**
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double> *in,
                               const unsigned int            *offsets,
                               double                        *out)
{
  const unsigned int n_chunks = n_entries/2;
  if (add_into)
    {
      for (unsigned int i=0; i<n_chunks; ++i)
        {
          __m128d u0 = in[2*i+0].data;
          __m128d u1 = in[2*i+1].data;
          __m128d res0 = _mm_unpacklo_pd (u0, u1);
          __m128d res1 = _mm_unpackhi_pd (u0, u1);
          _mm_storeu_pd(out+2*i+offsets[0], _mm_add_pd(_mm_loadu_pd(out+2*i+offsets[0]), res0));
          _mm_storeu_pd(out+2*i+offsets[1], _mm_add_pd(_mm_loadu_pd(out+2*i+offsets[1]), res1));
        }
      const unsigned int shift = n_chunks * 2;
      for (unsigned int i=shift; i<n_entries; ++i)
        for (unsigned int v=0; v<2; ++v)
          out[offsets[v]+i] += in[i][v];
    }
  else
    {
      for (unsigned int i=0; i<n_chunks; ++i)
        {
          __m128d u0 = in[2*i+0].data;
          __m128d u1 = in[2*i+1].data;
          __m128d res0 = _mm_unpacklo_pd (u0, u1);
          __m128d res1 = _mm_unpackhi_pd (u0, u1);
          _mm_storeu_pd(out+2*i+offsets[0], res0);
          _mm_storeu_pd(out+2*i+offsets[1], res1);
        }
      const unsigned int shift = n_chunks * 2;
      for (unsigned int i=shift; i<n_entries; ++i)
        for (unsigned int v=0; v<2; ++v)
          out[offsets[v]+i] = in[i][v];
    }
}



/**
 * Specialization for float and SSE2.
 */
template <>
class VectorizedArray<float>
{
public:
  /**
   * This gives the number of vectors collected in this class.
   */
  static const unsigned int n_array_elements = 4;

  /**
   * This function can be used to set all data fields to a given scalar.
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
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m128 data;

private:
  /**
   * Returns the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt () const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_ps(data);
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
    __m128 mask = _mm_set1_ps (-0.f);
    VectorizedArray res;
    res.data = _mm_andnot_ps(mask, data);
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
    res.data = _mm_max_ps (data, other.data);
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



/**
 * Specialization for float and SSE2.
 */
template <>
inline
void vectorized_load_and_transpose(const unsigned int      n_entries,
                                   const float            *in,
                                   const unsigned int     *offsets,
                                   VectorizedArray<float> *out)
{
  const unsigned int n_chunks = n_entries/4, remainder = n_entries%4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m128 u0 = _mm_loadu_ps(in+4*i+offsets[0]);
      __m128 u1 = _mm_loadu_ps(in+4*i+offsets[1]);
      __m128 u2 = _mm_loadu_ps(in+4*i+offsets[2]);
      __m128 u3 = _mm_loadu_ps(in+4*i+offsets[3]);
      __m128 v0 = _mm_shuffle_ps (u0, u1, 0x44);
      __m128 v1 = _mm_shuffle_ps (u0, u1, 0xee);
      __m128 v2 = _mm_shuffle_ps (u2, u3, 0x44);
      __m128 v3 = _mm_shuffle_ps (u2, u3, 0xee);
      out[4*i+0].data = _mm_shuffle_ps (v0, v2, 0x88);
      out[4*i+1].data = _mm_shuffle_ps (v0, v2, 0xdd);
      out[4*i+2].data = _mm_shuffle_ps (v1, v3, 0x88);
      out[4*i+3].data = _mm_shuffle_ps (v1, v3, 0xdd);
    }
  if (remainder > 0 && n_chunks > 0)
    {
      // simple re-load all data in the last slot
      const unsigned int final_pos = n_chunks*4-4+remainder;
      Assert(final_pos+4 == n_entries, ExcInternalError());
      __m128 u0 = _mm_loadu_ps(in+final_pos+offsets[0]);
      __m128 u1 = _mm_loadu_ps(in+final_pos+offsets[1]);
      __m128 u2 = _mm_loadu_ps(in+final_pos+offsets[2]);
      __m128 u3 = _mm_loadu_ps(in+final_pos+offsets[3]);
      __m128 v0 = _mm_shuffle_ps (u0, u1, 0x44);
      __m128 v1 = _mm_shuffle_ps (u0, u1, 0xee);
      __m128 v2 = _mm_shuffle_ps (u2, u3, 0x44);
      __m128 v3 = _mm_shuffle_ps (u2, u3, 0xee);
      out[final_pos+0].data = _mm_shuffle_ps (v0, v2, 0x88);
      out[final_pos+1].data = _mm_shuffle_ps (v0, v2, 0xdd);
      out[final_pos+2].data = _mm_shuffle_ps (v1, v3, 0x88);
      out[final_pos+3].data = _mm_shuffle_ps (v1, v3, 0xdd);
    }
  else if (remainder > 0)
    for (unsigned int i=0; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[i][v] = in[offsets[v]+i];
}



/**
 * Specialization for double and AVX.
 */
template <>
inline
void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float> *in,
                               const unsigned int           *offsets,
                               float                        *out)
{
  const unsigned int n_chunks = n_entries/4;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      __m128 u0 = in[4*i+0].data;
      __m128 u1 = in[4*i+1].data;
      __m128 u2 = in[4*i+2].data;
      __m128 u3 = in[4*i+3].data;
      __m128 t0 = _mm_shuffle_ps (u0, u1, 0x44);
      __m128 t1 = _mm_shuffle_ps (u0, u1, 0xee);
      __m128 t2 = _mm_shuffle_ps (u2, u3, 0x44);
      __m128 t3 = _mm_shuffle_ps (u2, u3, 0xee);
      u0 = _mm_shuffle_ps (t0, t2, 0x88);
      u1 = _mm_shuffle_ps (t0, t2, 0xdd);
      u2 = _mm_shuffle_ps (t1, t3, 0x88);
      u3 = _mm_shuffle_ps (t1, t3, 0xdd);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
        {
          u0 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[0]), u0);
          _mm_storeu_ps(out+4*i+offsets[0], u0);
          u1 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[1]), u1);
          _mm_storeu_ps(out+4*i+offsets[1], u1);
          u2 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[2]), u2);
          _mm_storeu_ps(out+4*i+offsets[2], u2);
          u3 = _mm_add_ps(_mm_loadu_ps(out+4*i+offsets[3]), u3);
          _mm_storeu_ps(out+4*i+offsets[3], u3);
        }
      else
        {
          _mm_storeu_ps(out+4*i+offsets[0], u0);
          _mm_storeu_ps(out+4*i+offsets[1], u1);
          _mm_storeu_ps(out+4*i+offsets[2], u2);
          _mm_storeu_ps(out+4*i+offsets[3], u3);
        }
    }
  const unsigned int shift = n_chunks * 4;
  if (add_into)
    for (unsigned int i=shift; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[offsets[v]+i] += in[i][v];
  else
    for (unsigned int i=shift; i<n_entries; ++i)
      for (unsigned int v=0; v<4; ++v)
        out[offsets[v]+i] = in[i][v];
}



#endif // if DEAL_II_COMPILER_VECTORIZATION_LEVEL > 0


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
    // read. This should save some instructions as compared to directly
    // setting the individual elements and also circumvents a compiler
    // optimization bug in gcc-4.6 with SSE2 (see also deal.II developers list
    // from April 2014, topic "matrix_free/step-48 Test").
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
   * Computes the tangent of a vectorized data field. The result is returned
   * as vectorized array in the form <tt>{tan(x[0]), tan(x[1]), ...,
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
   * Computes the exponential of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{exp(x[0]), exp(x[1]), ...,
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
   * Computes the square root of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{sqrt(x[0]), sqrt(x[1]),
   * ..., sqrt(x[n_array_elements-1])}</tt>.
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
   * Raises the given number @p x to the power @p p for a vectorized data
   * field. The result is returned as vectorized array in the form
   * <tt>{pow(x[0],p), pow(x[1],p), ..., pow(x[n_array_elements-1],p)}</tt>.
   *
   * @relates VectorizedArray
   */
  template <typename Number>
  inline
  ::dealii::VectorizedArray<Number>
  pow (const ::dealii::VectorizedArray<Number> &x,
       const Number p)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for (unsigned int i=0; i<dealii::VectorizedArray<Number>::n_array_elements; ++i)
      values[i] = std::pow(x[i], p);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
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
