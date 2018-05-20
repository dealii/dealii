// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

#ifndef dealii_vectorization_h
#define dealii_vectorization_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>

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

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 && !defined(__AVX__)
#  error \
    "Mismatch in vectorization capabilities: AVX was detected during configuration of deal.II and switched on, but it is apparently not available for the file you are trying to compile at the moment. Check compilation flags controlling the instruction set, such as -march=native."
#endif
#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 3 && !defined(__AVX512F__)
#  error \
    "Mismatch in vectorization capabilities: AVX-512F was detected during configuration of deal.II and switched on, but it is apparently not available for the file you are trying to compile at the moment. Check compilation flags controlling the instruction set, such as -march=native."
#endif

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 // AVX, AVX-512
#  include <immintrin.h>
#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL == 1 // SSE2
#  include <emmintrin.h>
#endif

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
  * The structs below are needed since VectorizedArray<T1> is a POD-type without a constructor and
  * can be a template argument for SymmetricTensor<...,T2> where T2 would equal VectorizedArray<T1>.
  * Internally, in previous versions of deal.II, SymmetricTensor<...,T2> would make use of the constructor
  * of T2 leading to a compile-time error. However simply adding a constructor for VectorizedArray<T1>
  * breaks the POD-idioms needed elsewhere. Calls to constructors of T2 subsequently got replaced by a
  * call to internal::NumberType<T2> which then determines the right function to use by template deduction.
  * A detailed discussion can be found at https://github.com/dealii/dealii/pull/3967 . Also see numbers.h
  * for other specializations.
  */
  template <typename T>
  struct NumberType<VectorizedArray<T>>
  {
    static const VectorizedArray<T>&
    value(const VectorizedArray<T>& t)
    {
      return t;
    }

    static VectorizedArray<T>
    value(const T& t)
    {
      VectorizedArray<T> tmp;
      tmp = t;
      return tmp;
    }
  };
} // namespace internal

// Enable the EnableIfScalar type trait for VectorizedArray<Number> such
// that it can be used as a Number type in Tensor<rank,dim,Number>, etc.

template <typename Number>
struct EnableIfScalar<VectorizedArray<Number>>
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
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const Number scalar)
  {
    data = scalar;
    return *this;
  }

  /**
   * Access operator (only valid with component 0)
   */
  DEAL_II_ALWAYS_INLINE
  Number& operator[](const unsigned int comp)
  {
    (void) comp;
    AssertIndexRange(comp, 1);
    return data;
  }

  /**
   * Constant access operator (only valid with component 0)
   */
  DEAL_II_ALWAYS_INLINE
  const Number& operator[](const unsigned int comp) const
  {
    (void) comp;
    AssertIndexRange(comp, 1);
    return data;
  }

  /**
   * Addition
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray<Number>& vec)
  {
    data += vec.data;
    return *this;
  }

  /**
   * Subtraction
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray<Number>& vec)
  {
    data -= vec.data;
    return *this;
  }

  /**
   * Multiplication
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray<Number>& vec)
  {
    data *= vec.data;
    return *this;
  }

  /**
   * Division
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray<Number>& vec)
  {
    data /= vec.data;
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by the amount of bytes
   * in the vectorized array, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const Number* ptr)
  {
    data = *ptr;
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * the amount of bytes in the vectorized array, as opposed to casting a
   * double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(Number* ptr) const
  {
    *ptr = data;
  }

  /**
   * Write the content of the calling class into memory in form of
   * @p n_array_elements to the given address using non-temporal stores that
   * bypass the processor's caches, using @p _mm_stream_pd store intrinsics on
   * supported CPUs. The destination of the store @p ptr must be aligned by
   * the amount of bytes in the vectorized array.
   *
   * This store operation can be faster than usual store operations in case
   * the store is streaming because it avoids the read-for-ownership transfer
   * typically invoked in standard stores. This approximately works as follows
   * (see the literature on computer architecture for details): When an
   * algorithm stores some results to a memory address, a processor typically
   * wants to move it into some of its caches as it expects the data to be
   * re-used again at some point. Since caches are organized in lines of sizes
   * either 64 byte or 128 byte but writes are usually smaller, a processor
   * must first load in the destination cache line upon a write because only
   * part of the cache line is overwritten initially. If a series of stores
   * write data in a chunk bigger than any of its caches could handle, the
   * data finally has to be moved out from the caches to main memory. But
   * since all addressed have first been read, this doubles the load on main
   * memory, which can incur a performance penalty. Furthermore, the
   * organization of caches in a multicore context also requires reading an
   * address before something can be written to cache to that address, see
   * e.g. the <a href="https://en.wikipedia.org/wiki/MESI_protocol">Wikipedia
   * article on the MESI protocol</a> for details. The instruction underlying
   * this function call signals to the processor that these two prerequisites
   * on a store are relaxed: Firstly, one expects the whole cache line to be
   * overwritten (that the memory subsystem then handles appropriately), so no
   * need to first read the "remainder" of the cache line. Secondly, the data
   * behind that particular memory will not be subject to cache coherency
   * protocol as it will be in main memory both when the same processor wants
   * to access it again as well as any other processors in a multicore
   * chip. Due to this particular setup, any subsequent access to the data
   * written by this function will need to query main memory, which is slower
   * than an access from a cache both latency-wise and throughput-wise. Thus,
   * this command should only be used for large stores that will collectively
   * not fit into caches, as performance will be degraded otherwise. For a
   * typical use case, see also <a
   * href="https://blogs.fau.de/hager/archives/2103">this blog article</a>.
   *
   * Note that streaming stores are only available in the specialized SSE/AVX
   * classes of VectorizedArray of type @p double or @p float, not in the
   * generic base class.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(Number* ptr) const
  {
    *ptr = data;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const Number* base_ptr, const unsigned int* offsets)
  {
    data = base_ptr[offsets[0]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, Number* base_ptr) const
  {
    base_ptr[offsets[0]] = data;
  }

  /**
   * Actual data field. Since this class represents a POD data type, it is
   * declared public.
   */
  Number data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = std::sqrt(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    VectorizedArray res;
    res.data = std::fabs(data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = std::max(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = std::min(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Create a vectorized array that sets all entries in the array to the given
 * scalar.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             make_vectorized_array(const Number& u)
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
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const Number*            in,
                              const unsigned int*      offsets,
                              VectorizedArray<Number>* out)
{
  for(unsigned int i = 0; i < n_entries; ++i)
    for(unsigned int v = 0; v < VectorizedArray<Number>::n_array_elements; ++v)
      out[i][v] = in[offsets[v] + i];
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
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<Number>* in,
                               const unsigned int*            offsets,
                               Number*                        out)
{
  if(add_into)
    for(unsigned int i = 0; i < n_entries; ++i)
      for(unsigned int v = 0; v < VectorizedArray<Number>::n_array_elements;
          ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for(unsigned int i = 0; i < n_entries; ++i)
      for(unsigned int v = 0; v < VectorizedArray<Number>::n_array_elements;
          ++v)
        out[offsets[v] + i] = in[i][v];
}

// for safety, also check that __AVX512F__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#if DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 3 && defined(__AVX512F__)

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
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const double x)
  {
    data = _mm512_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  double& operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<double*>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double& operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<const double*>(&data) + comp);
  }

  /**
   * Addition.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray& vec)
  {
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#  else
    data = _mm512_add_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Subtraction.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#  else
    data = _mm512_sub_pd(data, vec.data);
#  endif
    return *this;
  }
  /**
   * Multiplication.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#  else
    data = _mm512_mul_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Division.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#  else
    data = _mm512_div_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 64 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double* ptr)
  {
    data = _mm512_loadu_pd(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 64 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double* ptr) const
  {
    _mm512_storeu_pd(ptr, data);
  }

  /** @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 64 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double* ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 64 == 0,
           ExcMessage("Memory not aligned"));
    _mm512_stream_pd(ptr, data);
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double* base_ptr, const unsigned int* offsets)
  {
    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256  index_val = _mm256_loadu_ps((const float*) offsets);
    const __m256i index     = *((__m256i*) (&index_val));
    data                    = _mm512_i32gather_pd(index, base_ptr, 8);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, double* base_ptr) const
  {
    for(unsigned int i = 0; i < 8; ++i)
      for(unsigned int j = i + 1; j < 8; ++j)
        Assert(offsets[i] != offsets[j],
               ExcMessage("Result of scatter undefined if two offset elements"
                          " point to the same position"));

    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256  index_val = _mm256_loadu_ps((const float*) offsets);
    const __m256i index     = *((__m256i*) (&index_val));
    _mm512_i32scatter_pd(base_ptr, index, data, 8);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m512d data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm512_sqrt_pd(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +. Since there is no andnot for AVX512, we interpret the data
    // as 64 bit integers and do the andnot on those types (note that andnot
    // is a bitwise operation so the data type does not matter)
    __m512d         mask = _mm512_set1_pd(-0.);
    VectorizedArray res;
    res.data = (__m512d) _mm512_andnot_epi64((__m512i) mask, (__m512i) data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm512_max_pd(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_pd(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Specialization for double and AVX-512.
 */
template <>
inline void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const double*            in,
                              const unsigned int*      offsets,
                              VectorizedArray<double>* out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int outer = 0; outer < 8; outer += 4)
    {
      const double* in0 = in + offsets[0 + outer];
      const double* in1 = in + offsets[1 + outer];
      const double* in2 = in + offsets[2 + outer];
      const double* in3 = in + offsets[3 + outer];

      for(unsigned int i = 0; i < n_chunks; ++i)
        {
          __m256d u0 = _mm256_loadu_pd(in0 + 4 * i);
          __m256d u1 = _mm256_loadu_pd(in1 + 4 * i);
          __m256d u2 = _mm256_loadu_pd(in2 + 4 * i);
          __m256d u3 = _mm256_loadu_pd(in3 + 4 * i);
          __m256d t0 = _mm256_permute2f128_pd(u0, u2, 0x20);
          __m256d t1 = _mm256_permute2f128_pd(u1, u3, 0x20);
          __m256d t2 = _mm256_permute2f128_pd(u0, u2, 0x31);
          __m256d t3 = _mm256_permute2f128_pd(u1, u3, 0x31);
          *(__m256d*) ((double*) (&out[4 * i + 0].data) + outer)
            = _mm256_unpacklo_pd(t0, t1);
          *(__m256d*) ((double*) (&out[4 * i + 1].data) + outer)
            = _mm256_unpackhi_pd(t0, t1);
          *(__m256d*) ((double*) (&out[4 * i + 2].data) + outer)
            = _mm256_unpacklo_pd(t2, t3);
          *(__m256d*) ((double*) (&out[4 * i + 3].data) + outer)
            = _mm256_unpackhi_pd(t2, t3);
        }
      for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
        for(unsigned int v = 0; v < 4; ++v)
          out[i][outer + v] = in[offsets[v + outer] + i];
    }
}

/**
 * Specialization for double and AVX-512.
 */
template <>
inline void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double>* in,
                               const unsigned int*            offsets,
                               double*                        out)
{
  const unsigned int n_chunks = n_entries / 4;
  // do not do full transpose because the code is too long and will most
  // likely not pay off. rather do the transposition on the vectorized array
  // on size smaller, mm256d
  for(unsigned int outer = 0; outer < 8; outer += 4)
    {
      double* out0 = out + offsets[0 + outer];
      double* out1 = out + offsets[1 + outer];
      double* out2 = out + offsets[2 + outer];
      double* out3 = out + offsets[3 + outer];
      for(unsigned int i = 0; i < n_chunks; ++i)
        {
          __m256d u0
            = *(const __m256d*) ((const double*) (&in[4 * i + 0].data) + outer);
          __m256d u1
            = *(const __m256d*) ((const double*) (&in[4 * i + 1].data) + outer);
          __m256d u2
            = *(const __m256d*) ((const double*) (&in[4 * i + 2].data) + outer);
          __m256d u3
            = *(const __m256d*) ((const double*) (&in[4 * i + 3].data) + outer);
          __m256d t0   = _mm256_permute2f128_pd(u0, u2, 0x20);
          __m256d t1   = _mm256_permute2f128_pd(u1, u3, 0x20);
          __m256d t2   = _mm256_permute2f128_pd(u0, u2, 0x31);
          __m256d t3   = _mm256_permute2f128_pd(u1, u3, 0x31);
          __m256d res0 = _mm256_unpacklo_pd(t0, t1);
          __m256d res1 = _mm256_unpackhi_pd(t0, t1);
          __m256d res2 = _mm256_unpacklo_pd(t2, t3);
          __m256d res3 = _mm256_unpackhi_pd(t2, t3);

          // Cannot use the same store instructions in both paths of the 'if'
          // because the compiler cannot know that there is no aliasing between
          // pointers
          if(add_into)
            {
              res0 = _mm256_add_pd(_mm256_loadu_pd(out0 + 4 * i), res0);
              _mm256_storeu_pd(out0 + 4 * i, res0);
              res1 = _mm256_add_pd(_mm256_loadu_pd(out1 + 4 * i), res1);
              _mm256_storeu_pd(out1 + 4 * i, res1);
              res2 = _mm256_add_pd(_mm256_loadu_pd(out2 + 4 * i), res2);
              _mm256_storeu_pd(out2 + 4 * i, res2);
              res3 = _mm256_add_pd(_mm256_loadu_pd(out3 + 4 * i), res3);
              _mm256_storeu_pd(out3 + 4 * i, res3);
            }
          else
            {
              _mm256_storeu_pd(out0 + 4 * i, res0);
              _mm256_storeu_pd(out1 + 4 * i, res1);
              _mm256_storeu_pd(out2 + 4 * i, res2);
              _mm256_storeu_pd(out3 + 4 * i, res3);
            }
        }
      if(add_into)
        for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
          for(unsigned int v = 0; v < 4; ++v)
            out[offsets[v + outer] + i] += in[i][v + outer];
      else
        for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
          for(unsigned int v = 0; v < 4; ++v)
            out[offsets[v + outer] + i] = in[i][v + outer];
    }
}

/**
 * Specialization for float and AVX512.
 */
template <>
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
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const float x)
  {
    data = _mm512_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  float& operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 16);
    return *(reinterpret_cast<float*>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float& operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 16);
    return *(reinterpret_cast<const float*>(&data) + comp);
  }

  /**
   * Addition.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray& vec)
  {
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#  else
    data = _mm512_add_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Subtraction.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#  else
    data = _mm512_sub_ps(data, vec.data);
#  endif
    return *this;
  }
  /**
   * Multiplication.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#  else
    data = _mm512_mul_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Division.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#  else
    data = _mm512_div_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 64 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float* ptr)
  {
    data = _mm512_loadu_ps(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 64 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float* ptr) const
  {
    _mm512_storeu_ps(ptr, data);
  }

  /** @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 64 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float* ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 64 == 0,
           ExcMessage("Memory not aligned"));
    _mm512_stream_ps(ptr, data);
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float* base_ptr, const unsigned int* offsets)
  {
    // unfortunately, there does not appear to be a 512 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m512  index_val = _mm512_loadu_ps((const float*) offsets);
    const __m512i index     = *((__m512i*) (&index_val));
    data                    = _mm512_i32gather_ps(index, base_ptr, 4);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, float* base_ptr) const
  {
    for(unsigned int i = 0; i < 16; ++i)
      for(unsigned int j = i + 1; j < 16; ++j)
        Assert(offsets[i] != offsets[j],
               ExcMessage("Result of scatter undefined if two offset elements"
                          " point to the same position"));

    // unfortunately, there does not appear to be a 512 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m512  index_val = _mm512_loadu_ps((const float*) offsets);
    const __m512i index     = *((__m512i*) (&index_val));
    _mm512_i32scatter_ps(base_ptr, index, data, 4);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m512 data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm512_sqrt_ps(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +. Since there is no andnot for AVX512, we interpret the data
    // as 32 bit integers and do the andnot on those types (note that andnot
    // is a bitwise operation so the data type does not matter)
    __m512          mask = _mm512_set1_ps(-0.f);
    VectorizedArray res;
    res.data = (__m512) _mm512_andnot_epi32((__m512i) mask, (__m512i) data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm512_max_ps(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_ps(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Specialization for float and AVX-512.
 */
template <>
inline void
vectorized_load_and_transpose(const unsigned int      n_entries,
                              const float*            in,
                              const unsigned int*     offsets,
                              VectorizedArray<float>* out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int outer = 0; outer < 16; outer += 8)
    {
      for(unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128 u0 = _mm_loadu_ps(in + 4 * i + offsets[0 + outer]);
          __m128 u1 = _mm_loadu_ps(in + 4 * i + offsets[1 + outer]);
          __m128 u2 = _mm_loadu_ps(in + 4 * i + offsets[2 + outer]);
          __m128 u3 = _mm_loadu_ps(in + 4 * i + offsets[3 + outer]);
          __m128 u4 = _mm_loadu_ps(in + 4 * i + offsets[4 + outer]);
          __m128 u5 = _mm_loadu_ps(in + 4 * i + offsets[5 + outer]);
          __m128 u6 = _mm_loadu_ps(in + 4 * i + offsets[6 + outer]);
          __m128 u7 = _mm_loadu_ps(in + 4 * i + offsets[7 + outer]);
          // To avoid warnings about uninitialized variables, need to initialize
          // one variable with zero before using it.
          __m256 t0, t1, t2, t3 = _mm256_set1_ps(0.F);
          t0        = _mm256_insertf128_ps(t3, u0, 0);
          t0        = _mm256_insertf128_ps(t0, u4, 1);
          t1        = _mm256_insertf128_ps(t3, u1, 0);
          t1        = _mm256_insertf128_ps(t1, u5, 1);
          t2        = _mm256_insertf128_ps(t3, u2, 0);
          t2        = _mm256_insertf128_ps(t2, u6, 1);
          t3        = _mm256_insertf128_ps(t3, u3, 0);
          t3        = _mm256_insertf128_ps(t3, u7, 1);
          __m256 v0 = _mm256_shuffle_ps(t0, t1, 0x44);
          __m256 v1 = _mm256_shuffle_ps(t0, t1, 0xee);
          __m256 v2 = _mm256_shuffle_ps(t2, t3, 0x44);
          __m256 v3 = _mm256_shuffle_ps(t2, t3, 0xee);
          *(__m256*) ((float*) (&out[4 * i + 0].data) + outer)
            = _mm256_shuffle_ps(v0, v2, 0x88);
          *(__m256*) ((float*) (&out[4 * i + 1].data) + outer)
            = _mm256_shuffle_ps(v0, v2, 0xdd);
          *(__m256*) ((float*) (&out[4 * i + 2].data) + outer)
            = _mm256_shuffle_ps(v1, v3, 0x88);
          *(__m256*) ((float*) (&out[4 * i + 3].data) + outer)
            = _mm256_shuffle_ps(v1, v3, 0xdd);
        }
      for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
        for(unsigned int v = 0; v < 8; ++v)
          out[i][v + outer] = in[offsets[v + outer] + i];
    }
}

/**
 * Specialization for float and AVX-512.
 */
template <>
inline void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float>* in,
                               const unsigned int*           offsets,
                               float*                        out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int outer = 0; outer < 16; outer += 8)
    {
      for(unsigned int i = 0; i < n_chunks; ++i)
        {
          __m256 u0
            = *(const __m256*) ((const float*) (&in[4 * i + 0].data) + outer);
          __m256 u1
            = *(const __m256*) ((const float*) (&in[4 * i + 1].data) + outer);
          __m256 u2
            = *(const __m256*) ((const float*) (&in[4 * i + 2].data) + outer);
          __m256 u3
            = *(const __m256*) ((const float*) (&in[4 * i + 3].data) + outer);
          __m256 t0   = _mm256_shuffle_ps(u0, u1, 0x44);
          __m256 t1   = _mm256_shuffle_ps(u0, u1, 0xee);
          __m256 t2   = _mm256_shuffle_ps(u2, u3, 0x44);
          __m256 t3   = _mm256_shuffle_ps(u2, u3, 0xee);
          u0          = _mm256_shuffle_ps(t0, t2, 0x88);
          u1          = _mm256_shuffle_ps(t0, t2, 0xdd);
          u2          = _mm256_shuffle_ps(t1, t3, 0x88);
          u3          = _mm256_shuffle_ps(t1, t3, 0xdd);
          __m128 res0 = _mm256_extractf128_ps(u0, 0);
          __m128 res4 = _mm256_extractf128_ps(u0, 1);
          __m128 res1 = _mm256_extractf128_ps(u1, 0);
          __m128 res5 = _mm256_extractf128_ps(u1, 1);
          __m128 res2 = _mm256_extractf128_ps(u2, 0);
          __m128 res6 = _mm256_extractf128_ps(u2, 1);
          __m128 res3 = _mm256_extractf128_ps(u3, 0);
          __m128 res7 = _mm256_extractf128_ps(u3, 1);

          // Cannot use the same store instructions in both paths of the 'if'
          // because the compiler cannot know that there is no aliasing between
          // pointers
          if(add_into)
            {
              res0 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[0 + outer]),
                                res0);
              _mm_storeu_ps(out + 4 * i + offsets[0 + outer], res0);
              res1 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[1 + outer]),
                                res1);
              _mm_storeu_ps(out + 4 * i + offsets[1 + outer], res1);
              res2 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[2 + outer]),
                                res2);
              _mm_storeu_ps(out + 4 * i + offsets[2 + outer], res2);
              res3 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[3 + outer]),
                                res3);
              _mm_storeu_ps(out + 4 * i + offsets[3 + outer], res3);
              res4 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[4 + outer]),
                                res4);
              _mm_storeu_ps(out + 4 * i + offsets[4 + outer], res4);
              res5 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[5 + outer]),
                                res5);
              _mm_storeu_ps(out + 4 * i + offsets[5 + outer], res5);
              res6 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[6 + outer]),
                                res6);
              _mm_storeu_ps(out + 4 * i + offsets[6 + outer], res6);
              res7 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[7 + outer]),
                                res7);
              _mm_storeu_ps(out + 4 * i + offsets[7 + outer], res7);
            }
          else
            {
              _mm_storeu_ps(out + 4 * i + offsets[0 + outer], res0);
              _mm_storeu_ps(out + 4 * i + offsets[1 + outer], res1);
              _mm_storeu_ps(out + 4 * i + offsets[2 + outer], res2);
              _mm_storeu_ps(out + 4 * i + offsets[3 + outer], res3);
              _mm_storeu_ps(out + 4 * i + offsets[4 + outer], res4);
              _mm_storeu_ps(out + 4 * i + offsets[5 + outer], res5);
              _mm_storeu_ps(out + 4 * i + offsets[6 + outer], res6);
              _mm_storeu_ps(out + 4 * i + offsets[7 + outer], res7);
            }
        }
      if(add_into)
        for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
          for(unsigned int v = 0; v < 8; ++v)
            out[offsets[v + outer] + i] += in[i][v + outer];
      else
        for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
          for(unsigned int v = 0; v < 8; ++v)
            out[offsets[v + outer] + i] = in[i][v + outer];
    }
}

#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 2 && defined(__AVX__)

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
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const double x)
  {
    data = _mm256_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  double& operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<double*>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double& operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const double*>(&data) + comp);
  }

  /**
   * Addition.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray& vec)
  {
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#  else
    data = _mm256_add_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Subtraction.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#  else
    data = _mm256_sub_pd(data, vec.data);
#  endif
    return *this;
  }
  /**
   * Multiplication.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#  else
    data = _mm256_mul_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Division.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#  else
    data = _mm256_div_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 32 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double* ptr)
  {
    data = _mm256_loadu_pd(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 32 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double* ptr) const
  {
    _mm256_storeu_pd(ptr, data);
  }

  /** @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 32 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double* ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 32 == 0,
           ExcMessage("Memory not aligned"));
    _mm256_stream_pd(ptr, data);
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double* base_ptr, const unsigned int* offsets)
  {
#  ifdef __AVX2__
    // unfortunately, there does not appear to be a 128 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m128  index_val = _mm_loadu_ps((const float*) offsets);
    const __m128i index     = *((__m128i*) (&index_val));
    data                    = _mm256_i32gather_pd(base_ptr, index, 8);
#  else
    for(unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<double*>(&data) + i) = base_ptr[offsets[i]];
#  endif
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, double* base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for(unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double*>(&data) + i);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m256d data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_pd(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m256d         mask = _mm256_set1_pd(-0.);
    VectorizedArray res;
    res.data = _mm256_andnot_pd(mask, data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm256_max_pd(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_pd(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Specialization for double and AVX.
 */
template <>
inline void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const double*            in,
                              const unsigned int*      offsets,
                              VectorizedArray<double>* out)
{
  const unsigned int n_chunks = n_entries / 4;
  const double*      in0      = in + offsets[0];
  const double*      in1      = in + offsets[1];
  const double*      in2      = in + offsets[2];
  const double*      in3      = in + offsets[3];

  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256d u0          = _mm256_loadu_pd(in0 + 4 * i);
      __m256d u1          = _mm256_loadu_pd(in1 + 4 * i);
      __m256d u2          = _mm256_loadu_pd(in2 + 4 * i);
      __m256d u3          = _mm256_loadu_pd(in3 + 4 * i);
      __m256d t0          = _mm256_permute2f128_pd(u0, u2, 0x20);
      __m256d t1          = _mm256_permute2f128_pd(u1, u3, 0x20);
      __m256d t2          = _mm256_permute2f128_pd(u0, u2, 0x31);
      __m256d t3          = _mm256_permute2f128_pd(u1, u3, 0x31);
      out[4 * i + 0].data = _mm256_unpacklo_pd(t0, t1);
      out[4 * i + 1].data = _mm256_unpackhi_pd(t0, t1);
      out[4 * i + 2].data = _mm256_unpacklo_pd(t2, t3);
      out[4 * i + 3].data = _mm256_unpackhi_pd(t2, t3);
    }
  for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for(unsigned int v = 0; v < 4; ++v)
      out[i][v] = in[offsets[v] + i];
}

/**
 * Specialization for double and AVX.
 */
template <>
inline void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double>* in,
                               const unsigned int*            offsets,
                               double*                        out)
{
  const unsigned int n_chunks = n_entries / 4;
  double*            out0     = out + offsets[0];
  double*            out1     = out + offsets[1];
  double*            out2     = out + offsets[2];
  double*            out3     = out + offsets[3];
  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256d u0   = in[4 * i + 0].data;
      __m256d u1   = in[4 * i + 1].data;
      __m256d u2   = in[4 * i + 2].data;
      __m256d u3   = in[4 * i + 3].data;
      __m256d t0   = _mm256_permute2f128_pd(u0, u2, 0x20);
      __m256d t1   = _mm256_permute2f128_pd(u1, u3, 0x20);
      __m256d t2   = _mm256_permute2f128_pd(u0, u2, 0x31);
      __m256d t3   = _mm256_permute2f128_pd(u1, u3, 0x31);
      __m256d res0 = _mm256_unpacklo_pd(t0, t1);
      __m256d res1 = _mm256_unpackhi_pd(t0, t1);
      __m256d res2 = _mm256_unpacklo_pd(t2, t3);
      __m256d res3 = _mm256_unpackhi_pd(t2, t3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if(add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out0 + 4 * i), res0);
          _mm256_storeu_pd(out0 + 4 * i, res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out1 + 4 * i), res1);
          _mm256_storeu_pd(out1 + 4 * i, res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out2 + 4 * i), res2);
          _mm256_storeu_pd(out2 + 4 * i, res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out3 + 4 * i), res3);
          _mm256_storeu_pd(out3 + 4 * i, res3);
        }
      else
        {
          _mm256_storeu_pd(out0 + 4 * i, res0);
          _mm256_storeu_pd(out1 + 4 * i, res1);
          _mm256_storeu_pd(out2 + 4 * i, res2);
          _mm256_storeu_pd(out3 + 4 * i, res3);
        }
    }
  if(add_into)
    for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for(unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for(unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] = in[i][v];
}

/**
 * Specialization for float and AVX.
 */
template <>
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
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const float x)
  {
    data = _mm256_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  float& operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<float*>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float& operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<const float*>(&data) + comp);
  }

  /**
   * Addition.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray& vec)
  {
    // if the compiler supports vector arithmetics, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#  else
    data = _mm256_add_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Subtraction.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#  else
    data = _mm256_sub_ps(data, vec.data);
#  endif
    return *this;
  }
  /**
   * Multiplication.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#  else
    data = _mm256_mul_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Division.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#  else
    data = _mm256_div_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 32 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float* ptr)
  {
    data = _mm256_loadu_ps(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 32 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float* ptr) const
  {
    _mm256_storeu_ps(ptr, data);
  }

  /** @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 32 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float* ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 32 == 0,
           ExcMessage("Memory not aligned"));
    _mm256_stream_ps(ptr, data);
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float* base_ptr, const unsigned int* offsets)
  {
#  ifdef __AVX2__
    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256  index_val = _mm256_loadu_ps((const float*) offsets);
    const __m256i index     = *((__m256i*) (&index_val));
    data                    = _mm256_i32gather_ps(base_ptr, index, 4);
#  else
    for(unsigned int i = 0; i < 8; ++i)
      *(reinterpret_cast<float*>(&data) + i) = base_ptr[offsets[i]];
#  endif
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, float* base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for(unsigned int i = 0; i < 8; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float*>(&data) + i);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m256 data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm256_sqrt_ps(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m256          mask = _mm256_set1_ps(-0.f);
    VectorizedArray res;
    res.data = _mm256_andnot_ps(mask, data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm256_max_ps(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_ps(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Specialization for float and AVX.
 */
template <>
inline void
vectorized_load_and_transpose(const unsigned int      n_entries,
                              const float*            in,
                              const unsigned int*     offsets,
                              VectorizedArray<float>* out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0 = _mm_loadu_ps(in + 4 * i + offsets[0]);
      __m128 u1 = _mm_loadu_ps(in + 4 * i + offsets[1]);
      __m128 u2 = _mm_loadu_ps(in + 4 * i + offsets[2]);
      __m128 u3 = _mm_loadu_ps(in + 4 * i + offsets[3]);
      __m128 u4 = _mm_loadu_ps(in + 4 * i + offsets[4]);
      __m128 u5 = _mm_loadu_ps(in + 4 * i + offsets[5]);
      __m128 u6 = _mm_loadu_ps(in + 4 * i + offsets[6]);
      __m128 u7 = _mm_loadu_ps(in + 4 * i + offsets[7]);
      // To avoid warnings about uninitialized variables, need to initialize
      // one variable with zero before using it.
      __m256 t0, t1, t2, t3 = _mm256_set1_ps(0.F);
      t0                  = _mm256_insertf128_ps(t3, u0, 0);
      t0                  = _mm256_insertf128_ps(t0, u4, 1);
      t1                  = _mm256_insertf128_ps(t3, u1, 0);
      t1                  = _mm256_insertf128_ps(t1, u5, 1);
      t2                  = _mm256_insertf128_ps(t3, u2, 0);
      t2                  = _mm256_insertf128_ps(t2, u6, 1);
      t3                  = _mm256_insertf128_ps(t3, u3, 0);
      t3                  = _mm256_insertf128_ps(t3, u7, 1);
      __m256 v0           = _mm256_shuffle_ps(t0, t1, 0x44);
      __m256 v1           = _mm256_shuffle_ps(t0, t1, 0xee);
      __m256 v2           = _mm256_shuffle_ps(t2, t3, 0x44);
      __m256 v3           = _mm256_shuffle_ps(t2, t3, 0xee);
      out[4 * i + 0].data = _mm256_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm256_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm256_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm256_shuffle_ps(v1, v3, 0xdd);
    }
  for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for(unsigned int v = 0; v < 8; ++v)
      out[i][v] = in[offsets[v] + i];
}

/**
 * Specialization for float and AVX.
 */
template <>
inline void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float>* in,
                               const unsigned int*           offsets,
                               float*                        out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256 u0   = in[4 * i + 0].data;
      __m256 u1   = in[4 * i + 1].data;
      __m256 u2   = in[4 * i + 2].data;
      __m256 u3   = in[4 * i + 3].data;
      __m256 t0   = _mm256_shuffle_ps(u0, u1, 0x44);
      __m256 t1   = _mm256_shuffle_ps(u0, u1, 0xee);
      __m256 t2   = _mm256_shuffle_ps(u2, u3, 0x44);
      __m256 t3   = _mm256_shuffle_ps(u2, u3, 0xee);
      u0          = _mm256_shuffle_ps(t0, t2, 0x88);
      u1          = _mm256_shuffle_ps(t0, t2, 0xdd);
      u2          = _mm256_shuffle_ps(t1, t3, 0x88);
      u3          = _mm256_shuffle_ps(t1, t3, 0xdd);
      __m128 res0 = _mm256_extractf128_ps(u0, 0);
      __m128 res4 = _mm256_extractf128_ps(u0, 1);
      __m128 res1 = _mm256_extractf128_ps(u1, 0);
      __m128 res5 = _mm256_extractf128_ps(u1, 1);
      __m128 res2 = _mm256_extractf128_ps(u2, 0);
      __m128 res6 = _mm256_extractf128_ps(u2, 1);
      __m128 res3 = _mm256_extractf128_ps(u3, 0);
      __m128 res7 = _mm256_extractf128_ps(u3, 1);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if(add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[0]), res0);
          _mm_storeu_ps(out + 4 * i + offsets[0], res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[1]), res1);
          _mm_storeu_ps(out + 4 * i + offsets[1], res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[2]), res2);
          _mm_storeu_ps(out + 4 * i + offsets[2], res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[3]), res3);
          _mm_storeu_ps(out + 4 * i + offsets[3], res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[4]), res4);
          _mm_storeu_ps(out + 4 * i + offsets[4], res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[5]), res5);
          _mm_storeu_ps(out + 4 * i + offsets[5], res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[6]), res6);
          _mm_storeu_ps(out + 4 * i + offsets[6], res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[7]), res7);
          _mm_storeu_ps(out + 4 * i + offsets[7], res7);
        }
      else
        {
          _mm_storeu_ps(out + 4 * i + offsets[0], res0);
          _mm_storeu_ps(out + 4 * i + offsets[1], res1);
          _mm_storeu_ps(out + 4 * i + offsets[2], res2);
          _mm_storeu_ps(out + 4 * i + offsets[3], res3);
          _mm_storeu_ps(out + 4 * i + offsets[4], res4);
          _mm_storeu_ps(out + 4 * i + offsets[5], res5);
          _mm_storeu_ps(out + 4 * i + offsets[6], res6);
          _mm_storeu_ps(out + 4 * i + offsets[7], res7);
        }
    }
  if(add_into)
    for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for(unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for(unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] = in[i][v];
}

// for safety, also check that __SSE2__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#elif DEAL_II_COMPILER_VECTORIZATION_LEVEL >= 1

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
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const double x)
  {
    data = _mm_set1_pd(x);
    return *this;
  }

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  double& operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<double*>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double& operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<const double*>(&data) + comp);
  }

  /**
   * Addition.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#  else
    data = _mm_add_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Subtraction.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#  else
    data = _mm_sub_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Multiplication.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#  else
    data = _mm_mul_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Division.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#  else
    data = _mm_div_pd(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double* ptr)
  {
    data = _mm_loadu_pd(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double* ptr) const
  {
    _mm_storeu_pd(ptr, data);
  }

  /** @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 16 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double* ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    _mm_stream_pd(ptr, data);
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double* base_ptr, const unsigned int* offsets)
  {
    for(unsigned int i = 0; i < 2; ++i)
      *(reinterpret_cast<double*>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, double* base_ptr) const
  {
    for(unsigned int i = 0; i < 2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double*>(&data) + i);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m128d data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_pd(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform
    // bitwise andnot with -0. This will leave all
    // value and exponent bits unchanged but force
    // the sign value to +.
    __m128d         mask = _mm_set1_pd(-0.);
    VectorizedArray res;
    res.data = _mm_andnot_pd(mask, data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm_max_pd(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm_min_pd(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Specialization for double and SSE2.
 */
template <>
inline void
vectorized_load_and_transpose(const unsigned int       n_entries,
                              const double*            in,
                              const unsigned int*      offsets,
                              VectorizedArray<double>* out)
{
  const unsigned int n_chunks = n_entries / 2;
  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128d u0          = _mm_loadu_pd(in + 2 * i + offsets[0]);
      __m128d u1          = _mm_loadu_pd(in + 2 * i + offsets[1]);
      out[2 * i + 0].data = _mm_unpacklo_pd(u0, u1);
      out[2 * i + 1].data = _mm_unpackhi_pd(u0, u1);
    }
  for(unsigned int i = 2 * n_chunks; i < n_entries; ++i)
    for(unsigned int v = 0; v < 2; ++v)
      out[i][v] = in[offsets[v] + i];
}

/**
 * Specialization for double and SSE2.
 */
template <>
inline void
vectorized_transpose_and_store(const bool                     add_into,
                               const unsigned int             n_entries,
                               const VectorizedArray<double>* in,
                               const unsigned int*            offsets,
                               double*                        out)
{
  const unsigned int n_chunks = n_entries / 2;
  if(add_into)
    {
      for(unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(
            out + 2 * i + offsets[0],
            _mm_add_pd(_mm_loadu_pd(out + 2 * i + offsets[0]), res0));
          _mm_storeu_pd(
            out + 2 * i + offsets[1],
            _mm_add_pd(_mm_loadu_pd(out + 2 * i + offsets[1]), res1));
        }
      for(unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for(unsigned int v = 0; v < 2; ++v)
          out[offsets[v] + i] += in[i][v];
    }
  else
    {
      for(unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out + 2 * i + offsets[0], res0);
          _mm_storeu_pd(out + 2 * i + offsets[1], res1);
        }
      for(unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for(unsigned int v = 0; v < 2; ++v)
          out[offsets[v] + i] = in[i][v];
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

  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator=(const float x)
  {
    data = _mm_set1_ps(x);
    return *this;
  }

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  float& operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<float*>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float& operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const float*>(&data) + comp);
  }

  /**
   * Addition.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator+=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#  else
    data = _mm_add_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Subtraction.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator-=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#  else
    data = _mm_sub_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Multiplication.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator*=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#  else
    data = _mm_mul_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Division.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray&
  operator/=(const VectorizedArray& vec)
  {
#  ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#  else
    data = _mm_div_ps(data, vec.data);
#  endif
    return *this;
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float* ptr)
  {
    data = _mm_loadu_ps(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float* ptr) const
  {
    _mm_storeu_ps(ptr, data);
  }

  /** @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 16 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float* ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    _mm_stream_ps(ptr, data);
  }

  /**
   * Load @p n_array_elements from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float* base_ptr, const unsigned int* offsets)
  {
    for(unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<float*>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * n_array_elements to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int* offsets, float* base_ptr) const
  {
    for(unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float*>(&data) + i);
  }

  /**
   * Actual data field. Since this class represents a POD data type, it
   * remains public.
   */
  __m128 data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = _mm_sqrt_ps(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_abs() const
  {
    // to compute the absolute value, perform bitwise andnot with -0. This
    // will leave all value and exponent bits unchanged but force the sign
    // value to +.
    __m128          mask = _mm_set1_ps(-0.f);
    VectorizedArray res;
    res.data = _mm_andnot_ps(mask, data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm_max_ps(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray& other) const
  {
    VectorizedArray res;
    res.data = _mm_min_ps(data, other.data);
    return res;
  }

  /**
   * Make a few functions friends.
   */
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::sqrt(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::abs(const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::max(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
  template <typename Number2>
  friend VectorizedArray<Number2>
  std::min(const VectorizedArray<Number2>&, const VectorizedArray<Number2>&);
};

/**
 * Specialization for float and SSE2.
 */
template <>
inline void
vectorized_load_and_transpose(const unsigned int      n_entries,
                              const float*            in,
                              const unsigned int*     offsets,
                              VectorizedArray<float>* out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0           = _mm_loadu_ps(in + 4 * i + offsets[0]);
      __m128 u1           = _mm_loadu_ps(in + 4 * i + offsets[1]);
      __m128 u2           = _mm_loadu_ps(in + 4 * i + offsets[2]);
      __m128 u3           = _mm_loadu_ps(in + 4 * i + offsets[3]);
      __m128 v0           = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 v1           = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 v2           = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 v3           = _mm_shuffle_ps(u2, u3, 0xee);
      out[4 * i + 0].data = _mm_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm_shuffle_ps(v1, v3, 0xdd);
    }
  for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for(unsigned int v = 0; v < 4; ++v)
      out[i][v] = in[offsets[v] + i];
}

/**
 * Specialization for float and SSE2.
 */
template <>
inline void
vectorized_transpose_and_store(const bool                    add_into,
                               const unsigned int            n_entries,
                               const VectorizedArray<float>* in,
                               const unsigned int*           offsets,
                               float*                        out)
{
  const unsigned int n_chunks = n_entries / 4;
  for(unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0 = in[4 * i + 0].data;
      __m128 u1 = in[4 * i + 1].data;
      __m128 u2 = in[4 * i + 2].data;
      __m128 u3 = in[4 * i + 3].data;
      __m128 t0 = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 t1 = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 t2 = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 t3 = _mm_shuffle_ps(u2, u3, 0xee);
      u0        = _mm_shuffle_ps(t0, t2, 0x88);
      u1        = _mm_shuffle_ps(t0, t2, 0xdd);
      u2        = _mm_shuffle_ps(t1, t3, 0x88);
      u3        = _mm_shuffle_ps(t1, t3, 0xdd);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if(add_into)
        {
          u0 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[0]), u0);
          _mm_storeu_ps(out + 4 * i + offsets[0], u0);
          u1 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[1]), u1);
          _mm_storeu_ps(out + 4 * i + offsets[1], u1);
          u2 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[2]), u2);
          _mm_storeu_ps(out + 4 * i + offsets[2], u2);
          u3 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[3]), u3);
          _mm_storeu_ps(out + 4 * i + offsets[3], u3);
        }
      else
        {
          _mm_storeu_ps(out + 4 * i + offsets[0], u0);
          _mm_storeu_ps(out + 4 * i + offsets[1], u1);
          _mm_storeu_ps(out + 4 * i + offsets[2], u2);
          _mm_storeu_ps(out + 4 * i + offsets[3], u3);
        }
    }
  if(add_into)
    for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for(unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for(unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for(unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] = in[i][v];
}

#endif // if DEAL_II_COMPILER_VECTORIZATION_LEVEL > 0

/**
 * Relational operator == for VectorizedArray
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE bool
operator==(const VectorizedArray<Number>& lhs,
           const VectorizedArray<Number>& rhs)
{
  for(unsigned int i = 0; i < VectorizedArray<Number>::n_array_elements; ++i)
    if(lhs[i] != rhs[i])
      return false;

  return true;
}

/**
 * Addition of two vectorized arrays with operator +.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator+(const VectorizedArray<Number>& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp = u;
  return tmp += v;
}

/**
 * Subtraction of two vectorized arrays with operator -.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator-(const VectorizedArray<Number>& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp = u;
  return tmp -= v;
}

/**
 * Multiplication of two vectorized arrays with operator *.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator*(const VectorizedArray<Number>& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp = u;
  return tmp *= v;
}

/**
 * Division of two vectorized arrays with operator /.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator/(const VectorizedArray<Number>& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp = u;
  return tmp /= v;
}

/**
 * Addition of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator+(const Number& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp += v;
}

/**
 * Addition of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator+(const double& u, const VectorizedArray<float>& v)
{
  VectorizedArray<float> tmp;
  tmp = u;
  return tmp += v;
}

/**
 * Addition of a vectorized array and a scalar (expanded to a vectorized array
 * with @p n_array_elements equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator+(const VectorizedArray<Number>& v, const Number& u)
{
  return u + v;
}

/**
 * Addition of a vectorized array and a scalar (expanded to a vectorized array
 * with @p n_array_elements equal entries) in case the scalar is a double
 * (needed in order to be able to write simple code with constants that are
 * usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator+(const VectorizedArray<float>& v, const double& u)
{
  return u + v;
}

/**
 * Subtraction of a vectorized array from a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator-(const Number& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp -= v;
}

/**
 * Subtraction of a vectorized array from a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator-(const double& u, const VectorizedArray<float>& v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp -= v;
}

/**
 * Subtraction of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) from a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator-(const VectorizedArray<Number>& v, const Number& u)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return v - tmp;
}

/**
 * Subtraction of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) from a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator-(const VectorizedArray<float>& v, const double& u)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return v - tmp;
}

/**
 * Multiplication of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator*(const Number& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp *= v;
}

/**
 * Multiplication of a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator*(const double& u, const VectorizedArray<float>& v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp *= v;
}

/**
 * Multiplication of a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator*(const VectorizedArray<Number>& v, const Number& u)
{
  return u * v;
}

/**
 * Multiplication of a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator*(const VectorizedArray<float>& v, const double& u)
{
  return u * v;
}

/**
 * Quotient between a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator/(const Number& u, const VectorizedArray<Number>& v)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return tmp /= v;
}

/**
 * Quotient between a scalar (expanded to a vectorized array with @p
 * n_array_elements equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator/(const double& u, const VectorizedArray<float>& v)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return tmp /= v;
}

/**
 * Quotient between a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator/(const VectorizedArray<Number>& v, const Number& u)
{
  VectorizedArray<Number> tmp;
  tmp = u;
  return v / tmp;
}

/**
 * Quotient between a vectorized array and a scalar (expanded to a vectorized
 * array with @p n_array_elements equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float>
                             operator/(const VectorizedArray<float>& v, const double& u)
{
  VectorizedArray<float> tmp;
  tmp = float(u);
  return v / tmp;
}

/**
 * Unary operator + on a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator+(const VectorizedArray<Number>& u)
{
  return u;
}

/**
 * Unary operator - on a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number>
                             operator-(const VectorizedArray<Number>& u)
{
  // to get a negative sign, subtract the input from zero (could also
  // multiply by -1, but this one is slightly simpler)
  return VectorizedArray<Number>() - u;
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
   * Compute the sine of a vectorized data field. The result is returned as
   * vectorized array in the form <tt>{sin(x[0]), sin(x[1]), ...,
   * sin(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  sin(const ::dealii::VectorizedArray<Number>& x)
  {
    // put values in an array and later read in that array with an unaligned
    // read. This should save some instructions as compared to directly
    // setting the individual elements and also circumvents a compiler
    // optimization bug in gcc-4.6 with SSE2 (see also deal.II developers list
    // from April 2014, topic "matrix_free/step-48 Test").
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for(unsigned int i = 0;
        i < dealii::VectorizedArray<Number>::n_array_elements;
        ++i)
      values[i] = std::sin(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }

  /**
   * Compute the cosine of a vectorized data field. The result is returned as
   * vectorized array in the form <tt>{cos(x[0]), cos(x[1]), ...,
   * cos(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  cos(const ::dealii::VectorizedArray<Number>& x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for(unsigned int i = 0;
        i < dealii::VectorizedArray<Number>::n_array_elements;
        ++i)
      values[i] = std::cos(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }

  /**
   * Compute the tangent of a vectorized data field. The result is returned
   * as vectorized array in the form <tt>{tan(x[0]), tan(x[1]), ...,
   * tan(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  tan(const ::dealii::VectorizedArray<Number>& x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for(unsigned int i = 0;
        i < dealii::VectorizedArray<Number>::n_array_elements;
        ++i)
      values[i] = std::tan(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }

  /**
   * Compute the exponential of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{exp(x[0]), exp(x[1]), ...,
   * exp(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  exp(const ::dealii::VectorizedArray<Number>& x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for(unsigned int i = 0;
        i < dealii::VectorizedArray<Number>::n_array_elements;
        ++i)
      values[i] = std::exp(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }

  /**
   * Compute the natural logarithm of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{log(x[0]), log(x[1]), ...,
   * log(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  log(const ::dealii::VectorizedArray<Number>& x)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for(unsigned int i = 0;
        i < dealii::VectorizedArray<Number>::n_array_elements;
        ++i)
      values[i] = std::log(x[i]);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }

  /**
   * Compute the square root of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{sqrt(x[0]), sqrt(x[1]),
   * ..., sqrt(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  sqrt(const ::dealii::VectorizedArray<Number>& x)
  {
    return x.get_sqrt();
  }

  /**
   * Raises the given number @p x to the power @p p for a vectorized data
   * field. The result is returned as vectorized array in the form
   * <tt>{pow(x[0],p), pow(x[1],p), ..., pow(x[n_array_elements-1],p)}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  pow(const ::dealii::VectorizedArray<Number>& x, const Number p)
  {
    Number values[::dealii::VectorizedArray<Number>::n_array_elements];
    for(unsigned int i = 0;
        i < dealii::VectorizedArray<Number>::n_array_elements;
        ++i)
      values[i] = std::pow(x[i], p);
    ::dealii::VectorizedArray<Number> out;
    out.load(&values[0]);
    return out;
  }

  /**
   * Compute the absolute value (modulus) of a vectorized data field. The
   * result is returned as vectorized array in the form <tt>{abs(x[0]),
   * abs(x[1]), ..., abs(x[n_array_elements-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  abs(const ::dealii::VectorizedArray<Number>& x)
  {
    return x.get_abs();
  }

  /**
   * Compute the componentwise maximum of two vectorized data fields. The
   * result is returned as vectorized array in the form <tt>{max(x[0],y[0]),
   * max(x[1],y[1]), ...}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  max(const ::dealii::VectorizedArray<Number>& x,
      const ::dealii::VectorizedArray<Number>& y)
  {
    return x.get_max(y);
  }

  /**
   * Compute the componentwise minimum of two vectorized data fields. The
   * result is returned as vectorized array in the form <tt>{min(x[0],y[0]),
   * min(x[1],y[1]), ...}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number>
  inline ::dealii::VectorizedArray<Number>
  min(const ::dealii::VectorizedArray<Number>& x,
      const ::dealii::VectorizedArray<Number>& y)
  {
    return x.get_min(y);
  }

} // namespace std

#endif
