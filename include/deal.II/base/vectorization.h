// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_vectorization_h
#define dealii_vectorization_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/template_constraints.h>

#include <algorithm>
#include <array>
#include <cmath>

// Note:
// The flag DEAL_II_VECTORIZATION_WIDTH_IN_BITS is essentially constructed
// according to the following scheme (on x86-based architectures)
// #ifdef __AVX512F__
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 512
// #elif defined (__AVX__)
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 256
// #elif defined (__SSE2__)
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 128
// #else
// #define DEAL_II_VECTORIZATION_WIDTH_IN_BITS 0
// #endif
// In addition to checking the flags __AVX512F__, __AVX__ and __SSE2__, a CMake
// test, 'check_01_cpu_features.cmake', ensures that these feature are not only
// present in the compilation unit but also working properly.

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS > 0

// These error messages try to detect the case that deal.II was compiled with
// a wider instruction set extension as the current compilation unit, for
// example because deal.II was compiled with AVX, but a user project does not
// add -march=native or similar flags, making it fall to SSE2. This leads to
// very strange errors as the size of data structures differs between the
// compiled deal.II code sitting in libdeal_II.so and the user code if not
// detected.
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && !defined(__AVX__)
#    error \
      "Mismatch in vectorization capabilities: AVX was detected during configuration of deal.II and switched on, but it is apparently not available for the file you are trying to compile at the moment. Check compilation flags controlling the instruction set, such as -march=native."
#  endif
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && !defined(__AVX512F__)
#    error \
      "Mismatch in vectorization capabilities: AVX-512F was detected during configuration of deal.II and switched on, but it is apparently not available for the file you are trying to compile at the moment. Check compilation flags controlling the instruction set, such as -march=native."
#  endif

#  ifdef _MSC_VER
#    include <intrin.h>
#  elif defined(__ALTIVEC__)
#    include <altivec.h>

// altivec.h defines vector, pixel, bool, but we do not use them, so undefine
// them before they make trouble
#    undef vector
#    undef pixel
#    undef bool
#  elif defined(__ARM_NEON)
#    include <arm_neon.h>
#  elif defined(__x86_64__)
#    include <x86intrin.h>
#  endif

#endif


DEAL_II_NAMESPACE_OPEN


// Enable the EnableIfScalar type trait for VectorizedArray<Number> such
// that it can be used as a Number type in Tensor<rank,dim,Number>, etc.

template <typename Number, std::size_t width>
struct EnableIfScalar<VectorizedArray<Number, width>>
{
  using type = VectorizedArray<typename EnableIfScalar<Number>::type, width>;
};



/**
 * An iterator for VectorizedArray.
 */
template <typename T>
class VectorizedArrayIterator
{
public:
  /**
   * Constructor.
   *
   * @param data The actual VectorizedArray.
   * @param lane A pointer to the current lane.
   */
  constexpr VectorizedArrayIterator(T &data, const std::size_t lane)
    : data(&data)
    , lane(lane)
  {}

  /**
   * Compare for equality.
   */
  constexpr bool
  operator==(const VectorizedArrayIterator<T> &other) const
  {
    Assert(this->data == other.data,
           ExcMessage(
             "You are trying to compare iterators into different arrays."));
    return this->lane == other.lane;
  }

  /**
   * Compare for inequality.
   */
  constexpr bool
  operator!=(const VectorizedArrayIterator<T> &other) const
  {
    Assert(this->data == other.data,
           ExcMessage(
             "You are trying to compare iterators into different arrays."));
    return this->lane != other.lane;
  }

  /**
   * Dereferencing operator (const version): returns the value of the current
   * lane.
   */
  constexpr const typename T::value_type &
  operator*() const
  {
    AssertIndexRange(lane, T::size());
    return (*data)[lane];
  }


  /**
   * Dereferencing operator (non-@p const version): returns the value of the
   * current lane.
   */
  template <typename U = T>
  constexpr std::enable_if_t<!std::is_same_v<U, const U>,
                             typename T::value_type> &
  operator*()
  {
    AssertIndexRange(lane, T::size());
    return (*data)[lane];
  }

  /**
   * Prefix <tt>++</tt> operator: <tt>++iterator</tt>. This operator advances
   * the iterator to the next lane and returns a reference to
   * <tt>*this</tt>.
   */
  constexpr VectorizedArrayIterator<T> &
  operator++()
  {
    AssertIndexRange(lane + 1, T::size() + 1);
    ++lane;
    return *this;
  }

  /**
   * This operator advances the iterator by @p offset lanes and returns a
   * reference to <tt>*this</tt>.
   */
  constexpr VectorizedArrayIterator<T> &
  operator+=(const std::size_t offset)
  {
    AssertIndexRange(lane + offset, T::size() + 1);
    lane += offset;
    return *this;
  }

  /**
   * Prefix <tt>--</tt> operator: <tt>--iterator</tt>. This operator advances
   * the iterator to the previous lane and returns a reference to
   * <tt>*this</tt>.
   */
  constexpr VectorizedArrayIterator<T> &
  operator--()
  {
    Assert(
      lane > 0,
      ExcMessage(
        "You can't decrement an iterator that is already at the beginning of the range."));
    --lane;
    return *this;
  }

  /**
   * Create new iterator, which is shifted by @p offset.
   */
  constexpr VectorizedArrayIterator<T>
  operator+(const std::size_t &offset) const
  {
    AssertIndexRange(lane + offset, T::size() + 1);
    return VectorizedArrayIterator<T>(*data, lane + offset);
  }

  /**
   * Compute distance between this iterator and iterator @p other.
   */
  constexpr std::ptrdiff_t
  operator-(const VectorizedArrayIterator<T> &other) const
  {
    return static_cast<std::ptrdiff_t>(lane) -
           static_cast<std::ptrdiff_t>(other.lane);
  }

private:
  /**
   * Pointer to the actual VectorizedArray.
   */
  T *data;

  /**
   * Pointer to the current lane.
   */
  std::size_t lane;
};



/**
 * A base class for the various VectorizedArray template specializations,
 * containing common functionalities.
 *
 * @tparam VectorizedArrayType Type of the actual vectorized array this
 *   class is operating on. We are using the
 *   Curiously Recurring Template Pattern (see
 *   https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) in this
 *   class to avoid having to resort to `virtual` member functions. In
 *   other words, `VectorizedArrayType` is a class derived from the
 *   current class.
 */
template <typename VectorizedArrayType, std::size_t width>
class VectorizedArrayBase
{
public:
  /**
   * Default constructor.
   */
  constexpr VectorizedArrayBase() = default;

  /**
   * Construct an array with the given initializer list.
   *
   * The initializer list must have at most as many elements as the
   * vector length. Elements not listed in the initializer list are
   * zero-initialized.
   */
  template <typename U>
  constexpr VectorizedArrayBase(const std::initializer_list<U> &list)
  {
    const unsigned int n_initializers = list.size();
    Assert(n_initializers <= size(),
           ExcMessage("The initializer list must have at most "
                      "as many elements as the vector length."));

    // Copy what's in the list.
    std::copy_n(list.begin(), n_initializers, this->begin());

    // Then add zero padding where necessary.
    if (n_initializers < size())
      std::fill(this->begin() + n_initializers, this->end(), 0.0);
  }

  /**
   * Return the number of elements in the array.
   */
  static constexpr std::size_t
  size()
  {
    return width;
  }

  /**
   * @return An iterator pointing to the beginning of the underlying data.
   */
  constexpr VectorizedArrayIterator<VectorizedArrayType>
  begin()
  {
    return VectorizedArrayIterator<VectorizedArrayType>(
      static_cast<VectorizedArrayType &>(*this), 0);
  }

  /**
   * @return An iterator pointing to the beginning of the underlying data (`const`
   * version).
   */
  constexpr VectorizedArrayIterator<const VectorizedArrayType>
  begin() const
  {
    return VectorizedArrayIterator<const VectorizedArrayType>(
      static_cast<const VectorizedArrayType &>(*this), 0);
  }

  /**
   * @return An iterator pointing to the end of the underlying data.
   */
  constexpr VectorizedArrayIterator<VectorizedArrayType>
  end()
  {
    return VectorizedArrayIterator<VectorizedArrayType>(
      static_cast<VectorizedArrayType &>(*this), width);
  }

  /**
   * @return An iterator pointing to the end of the underlying data (`const`
   * version).
   */
  constexpr VectorizedArrayIterator<const VectorizedArrayType>
  end() const
  {
    return VectorizedArrayIterator<const VectorizedArrayType>(
      static_cast<const VectorizedArrayType &>(*this), width);
  }

  /**
   * A default implementation for computing the dot product between two
   * vectorized arrays. It first forms the lane-by-lane product between
   * the current object and the argument, and then the across-lanes
   * sum of the product. The function returns the kind of object that
   * the VectorizedArray::sum() function returns.
   *
   * This function is inherited by all derived classes and provides
   * the dot product to them unless they override it with their own
   * implementation (presumably using a more efficient approach).
   */
  auto
  dot_product(const VectorizedArrayType &v) const
  {
    VectorizedArrayType p = static_cast<const VectorizedArrayType &>(*this);
    p *= v;
    return p.sum();
  }
};



/**
 * This generic class defines a unified interface to a vectorized data type.
 * For general template arguments, this class simply corresponds to the
 * template argument. For example, VectorizedArray<long double> is nothing
 * else but a wrapper around <tt>long double</tt> with exactly one data field
 * of type <tt>long double</tt> and overloaded arithmetic operations. This
 * means that <tt>VectorizedArray<ComplicatedType></tt> has a similar layout
 * as ComplicatedType, provided that ComplicatedType defines basic arithmetic
 * operations. For floats and doubles, an array of numbers are packed together
 * with the goal to be processed in a single-instruction/multiple-data (SIMD)
 * fashion. In the SIMD context, the elements of such a short vector are often
 * called "lanes". The number of elements packed together, i.e., the number of
 * lanes, depends on the computer system and compiler flags that are used for
 * compilation of deal.II. The fundamental idea of these packed data types is
 * to use one single CPU instruction to perform arithmetic operations on the
 * whole array using the processor's vector (SIMD) units. Most computer
 * systems by 2010 standards will use an array of two doubles or four floats,
 * respectively (this corresponds to the SSE/SSE2 data sets) when compiling
 * deal.II on 64-bit operating systems. On later processors (such as Intel Sandy
 * Bridge and newer, or AMD Bulldozer processors and newer), four
 * doubles or eight floats are used when deal.II is configured using gcc with
 * `--with-cpu=native` or `--with-cpu=corei7-avx`. On compilations with
 * AVX-512 support (e.g., Intel Skylake Server from 2017), eight doubles
 * or sixteen floats are used.
 *
 * The behavior of this class is made similar to the basic data types double
 * and float. The definition of a vectorized array does not initialize the
 * data field but rather leaves it undefined, as is the case for double and
 * float. However, when calling something like `VectorizedArray<double> a =
 * VectorizedArray<double>()` or `VectorizedArray<double> a = 0.`, it sets all
 * numbers in this field to zero. This class is of standard layout type
 * according to the C++11 standard, which means that there is an equivalent C
 * representation and the class can e.g. be safely copied with std::memcpy.
 * (See also https://en.cppreference.com/w/cpp/named_req/StandardLayoutType.)
 * The standard layout is also necessary for ensuring correct alignment of
 * data with address boundaries when collected in a vector (i.e., when the
 * first element in a vector is properly aligned, all subsequent elements will
 * be correctly aligned, too).
 *
 * Note that for proper functioning of this class, certain data alignment
 * rules must be respected. This is because the computer expects the starting
 * address of a VectorizedArray<double> field at specific addresses in memory
 * (usually, the address of the vectorized array should be a multiple of the
 * length of the array in bytes). Otherwise, a segmentation fault or a severe
 * loss of performance might occur. When creating a single data field on the
 * stack like `VectorizedArray<double> a = 5.;`, the compiler will take care
 * of data alignment automatically. However, when allocating a long vector of
 * VectorizedArray<double> data, one needs to respect these rules. Use the
 * class AlignedVector or data containers based on AlignedVector (such as
 * Table) for this purpose. It is a class very similar to std::vector
 * otherwise but always makes sure that data is correctly aligned.
 *
 * The user can explicitly control the width of a particular instruction set
 * architecture (ISA) extension by specifying the number of lanes via the second
 * template parameter of this wrapper class. For example on Intel Skylake
 * Server, you have the following options for the data type double:
 *  - VectorizedArray<double, 1> // no vectorization (auto-optimization)
 *  - VectorizedArray<double, 2> // SSE2
 *  - VectorizedArray<double, 4> // AVX
 *  - VectorizedArray<double, 8> // AVX-512 (default)
 *
 * and for Intel Sandy Bridge, Haswell, Broadwell, AMD Bulldozer and Zen/Ryzen:
 *  - VectorizedArray<double, 1> // no vectorization (auto-optimization)
 *  - VectorizedArray<double, 2> // SSE2
 *  - VectorizedArray<double, 4> // AVX (default)
 *
 * and for processors with AltiVec support:
 *  - VectorizedArray<double, 1>
 *  - VectorizedArray<double, 2>
 *
 * For older x86 processors or in case no processor-specific compilation flags
 * were added (i.e., without `-D CMAKE_CXX_FLAGS=-march=native` or similar
 * flags):
 *  - VectorizedArray<double, 1> // no vectorization (auto-optimization)
 *  - VectorizedArray<double, 2> // SSE2
 *
 * Similar considerations also apply to the data type `float`.
 *
 * Wrongly selecting the width, e.g., width=3 or width=8 on a processor which
 * does not support AVX-512 leads to a failing `static_assert`, i.e., you cannot
 * compile code using a VectorizedArray class with a number of lanes that
 * is not supported by the platform you are on. In particular, all platforms
 * we are aware of only ever support widths that are powers of two.
 *
 * @tparam Number The underlying scalar data type.
 * @tparam width  Vector length. (Optional; if not set, the maximal width of the
 *                architecture is used.)
 */
template <typename Number, std::size_t width>
class VectorizedArray
  : public VectorizedArrayBase<VectorizedArray<Number, width>, 1>
{
public:
  /**
   * The scalar type of the array elements.
   */
  using value_type = Number;

  /**
   * A constexpr boolean indicating whether the VectorizedArray with the
   * given choice of template parameters @p Number and @p width is indeed
   * implemented. The generic implementation is only implemented for
   * @p width equal to one. For specializations of this class (which are
   * defined depending on the instruction sets available) the boolean is
   * set to true as well.
   */
  static constexpr bool is_implemented = (width == 1);

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const Number scalar)
  {
    static_assert(width == 1,
                  "You specified an illegal width that is not supported.");

    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<Number, width>, 1>(list)
  {
    static_assert(width == 1,
                  "You specified an illegal width that is not supported.");
  }

  /**
   * This function assigns a scalar to the current object.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const Number scalar) &
  {
    data = scalar;
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const Number scalar) && = delete;

  /**
   * Access operator (only valid with component 0 in the base class without
   * specialization).
   */
  DEAL_II_ALWAYS_INLINE
  Number &
  operator[](const unsigned int comp)
  {
    (void)comp;
    AssertIndexRange(comp, 1);
    return data;
  }

  /**
   * Constant access operator (only valid with component 0 in the base class
   * without specialization).
   */
  DEAL_II_ALWAYS_INLINE
  const Number &
  operator[](const unsigned int comp) const
  {
    (void)comp;
    AssertIndexRange(comp, 1);
    return data;
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data += vec.data;
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data -= vec.data;
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data *= vec.data;
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data /= vec.data;
    return *this;
  }

  /**
   * Load size() data items from memory into the calling class, starting at
   * the given address. The pointer `ptr` needs not be aligned by the amount
   * of bytes in the vectorized array, as opposed to casting a double address
   * to VectorizedArray<double>*.
   */
  template <typename OtherNumber>
  DEAL_II_ALWAYS_INLINE void
  load(const OtherNumber *ptr)
  {
    data = *ptr;
  }

  /**
   * Write the content of the calling class into memory in form of
   * size() data items to the given address. The pointer `ptr` needs not be
   * aligned by the amount of bytes in the vectorized array, as opposed to
   * casting a double address to VectorizedArray<double>*.
   */
  template <typename OtherNumber>
  DEAL_II_ALWAYS_INLINE void
  store(OtherNumber *ptr) const
  {
    *ptr = data;
  }

  /**
   * Write the content of the calling class into memory in form of
   * size() data items to the given address using non-temporal stores that
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
   * overwritten (meaning that the memory subsystem makes sure that
   * consecutive stores that together span a cache line are merged, and
   * appropriately handling the case where only part of a cache line is
   * written), so there is no need to first read the "remainder" of the cache
   * line. Secondly, the data behind that particular memory will not be
   * subject to cache coherency protocol as it will be in main memory both
   * when the same processor wants to access it again as well as any other
   * processors in a multicore chip. Due to this particular setup, any
   * subsequent access to the data written by this function will need to query
   * main memory, which is slower than an access from a cache both
   * latency-wise and throughput-wise. Thus, this command should only be used
   * for storing large arrays that will collectively not fit into caches, as
   * performance will be degraded otherwise. For a typical use case, see also
   * <a href="https://blogs.fau.de/hager/archives/2103">this blog article</a>.
   *
   * Note that streaming stores are only available in the specialized SSE/AVX
   * classes of VectorizedArray of type @p double or @p float, not in the
   * generic base class.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(Number *ptr) const
  {
    *ptr = data;
  }

  /**
   * Load size() data items from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const Number *base_ptr, const unsigned int *offsets)
  {
    data = base_ptr[offsets[0]];
  }

  /**
   * Write the content of the calling class into memory in form of
   * size() data items to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, Number *base_ptr) const
  {
    base_ptr[offsets[0]] = data;
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  DEAL_II_ALWAYS_INLINE
  Number
  sum() const
  {
    return data;
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
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
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = std::min(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * @name Packing and unpacking of a VectorizedArray
 * @{
 */

/**
 * Create a vectorized array that sets all entries in the array to the given
 * scalar, i.e., broadcasts the scalar to all array elements.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number,
          std::size_t width =
            internal::VectorizedArrayWidthSpecifier<Number>::max_width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             make_vectorized_array(const Number &u)
{
  VectorizedArray<Number, width> result = u;
  return result;
}



/**
 * Create a vectorized array of given type and broadcast the scalar value
 * to all array elements.
 *
 * @relatesalso VectorizedArray
 */
template <typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
make_vectorized_array(const typename VectorizedArrayType::value_type &u)
{
  static_assert(
    std::is_same_v<VectorizedArrayType,
                   VectorizedArray<typename VectorizedArrayType::value_type,
                                   VectorizedArrayType::size()>>,
    "VectorizedArrayType is not a VectorizedArray.");

  VectorizedArrayType result = u;
  return result;
}



/**
 * Load size() data items from memory into the VectorizedArray @p out,
 * starting at the given addresses and with given offset, each entry from the
 * offset providing one element of the vectorized array.
 *
 * This operation corresponds to the following code:
 * @code
 * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *   out.data[v] = ptrs[v][offset];
 * @endcode
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
gather(VectorizedArray<Number, width>    &out,
       const std::array<Number *, width> &ptrs,
       const unsigned int                 offset)
{
  for (unsigned int v = 0; v < width; ++v)
    out.data[v] = ptrs[v][offset];
}



/**
 * This method loads VectorizedArray::size() data streams from the
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
 *   for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *     out[i][v] = in[offsets[v]+i];
 * @endcode
 *
 * A more optimized version of this code will be used for supported types.
 *
 * This is the inverse operation to vectorized_transpose_and_store().
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int              n_entries,
                              const Number                   *in,
                              const unsigned int             *offsets,
                              VectorizedArray<Number, width> *out)
{
  for (unsigned int i = 0; i < n_entries; ++i)
    for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
      out[i][v] = in[offsets[v] + i];
}


/**
 * The same as above with the difference that an array of pointers are
 * passed in as input argument @p in.
 *
 * In analogy to the function above, one can consider that
 * `in+offset[v]` is precomputed and passed as input argument.
 *
 * However, this function can also be used if some function returns an array
 * of pointers and no assumption can be made that they belong to the same array,
 * i.e., they can have their origin in different memory allocations.
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int                 n_entries,
                              const std::array<Number *, width> &in,
                              VectorizedArray<Number, width>    *out)
{
  for (unsigned int i = 0; i < n_entries; ++i)
    for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
      out[i][v] = in[v][i];
}



/**
 * This method stores the vectorized arrays in transposed form into the given
 * output array @p out with the given offsets @p offsets. This operation
 * corresponds to a transformation of a struct-of-array (input) into an
 * array-of-struct (output). This method operates on plain array, so no checks
 * for valid data access are made. It is the user's responsibility to ensure
 * that the given arrays are valid according to the access layout below.
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
 *   for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *     out[offsets[v]+i] = in[i][v];
 * @endcode
 *
 * For <code>add_into == true</code>, the code implements the following
 * action:
 * @code
 * for (unsigned int i=0; i<n_entries; ++i)
 *   for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
 *     out[offsets[v]+i] += in[i][v];
 * @endcode
 *
 * A more optimized version of this code will be used for supported types.
 *
 * This is the inverse operation to vectorized_load_and_transpose().
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                            add_into,
                               const unsigned int                    n_entries,
                               const VectorizedArray<Number, width> *in,
                               const unsigned int                   *offsets,
                               Number                               *out)
{
  if (add_into)
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[offsets[v] + i] = in[i][v];
}


/**
 * The same as above with the difference that an array of pointers are
 * passed in as input argument @p out.
 *
 * In analogy to the function above, one can consider that
 * `out+offset[v]` is precomputed and passed as input argument.
 *
 * However, this function can also be used if some function returns an array
 * of pointers and no assumption can be made that they belong to the same array,
 * i.e., they can have their origin in different memory allocations.
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                            add_into,
                               const unsigned int                    n_entries,
                               const VectorizedArray<Number, width> *in,
                               std::array<Number *, width>          &out)
{
  if (add_into)
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 0; i < n_entries; ++i)
      for (unsigned int v = 0; v < VectorizedArray<Number, width>::size(); ++v)
        out[v][i] = in[i][v];
}


/** @} */

#ifndef DOXYGEN

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ARM_NEON)

/**
 * Specialization for double and ARM Neon.
 */
template <>
class VectorizedArray<double, 2>
  : public VectorizedArrayBase<VectorizedArray<double, 2>, 2>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = double;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<double, 2>, 2>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator=(const double x) &
  {
    data = vdupq_n_f64(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const double scalar) && = delete;

  /**
   * Access operator.
   */
  double &
  operator[](const unsigned int comp)
  {
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  const double &
  operator[](const unsigned int comp) const
  {
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data = vaddq_f64(data, vec.data);
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data = vsubq_f64(data, vec.data);
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data = vmulq_f64(data, vec.data);
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data = vdivq_f64(data, vec.data);
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  void
  load(const double *ptr)
  {
    data = vld1q_f64(ptr);
  }

  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    DEAL_II_OPENMP_SIMD_PRAGMA
    for (unsigned int i = 0; i < 2; ++i)
      data[i] = ptr[i];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  void
  store(double *ptr) const
  {
    vst1q_f64(ptr, data);
  }

  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    DEAL_II_OPENMP_SIMD_PRAGMA
    for (unsigned int i = 0; i < 2; ++i)
      ptr[i] = data[i];
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 16 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    vst1q_f64(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 2; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    for (unsigned int i = 0; i < 2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  double
  sum() const
  {
    return vaddvq_f64(data);
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  mutable float64x2_t data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = vsqrtq_f64(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs() const
  {
    VectorizedArray res;
    res.data = vabsq_f64(data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vmaxq_f64(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vminq_f64(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};

/**
 * Specialization for float and ARM Neon.
 */
template <>
class VectorizedArray<float, 4>
  : public VectorizedArrayBase<VectorizedArray<float, 4>, 4>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = float;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<float, 4>, 4>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  VectorizedArray &
  operator=(const float x) &
  {
    data = vdupq_n_f32(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const float scalar) && = delete;

  /**
   * Access operator.
   */
  value_type &
  operator[](const unsigned int comp)
  {
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  const value_type &
  operator[](const unsigned int comp) const
  {
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data = vaddq_f32(data, vec.data);
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data = vsubq_f32(data, vec.data);
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data = vmulq_f32(data, vec.data);
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data = vdivq_f32(data, vec.data);
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  void
  load(const float *ptr)
  {
    data = vld1q_f32(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  void
  store(float *ptr) const
  {
    vst1q_f32(ptr, data);
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 16 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    vst1q_f32(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  float
  sum() const
  {
    return vaddvq_f32(data);
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  mutable float32x4_t data;

private:
  /**
   * Return the square root of this field. Not for use in user code. Use
   * sqrt(x) instead.
   */
  VectorizedArray
  get_sqrt() const
  {
    VectorizedArray res;
    res.data = vsqrtq_f32(data);
    return res;
  }

  /**
   * Return the absolute value of this field. Not for use in user code. Use
   * abs(x) instead.
   */
  VectorizedArray
  get_abs() const
  {
    VectorizedArray res;
    res.data = vabsq_f32(data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vmaxq_f32(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vminq_f32(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};


#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)

/**
 * Specialization for double and SSE2.
 */
template <>
class VectorizedArray<double, 2>
  : public VectorizedArrayBase<VectorizedArray<double, 2>, 2>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = double;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<double, 2>, 2>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x) &
  {
    data = _mm_set1_pd(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const double scalar) && = delete;

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  double &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm_add_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm_sub_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm_mul_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm_div_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = _mm_loadu_pd(ptr);
  }

  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    DEAL_II_OPENMP_SIMD_PRAGMA
    for (unsigned int i = 0; i < 2; ++i)
      data[i] = ptr[i];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    _mm_storeu_pd(ptr, data);
  }

  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    DEAL_II_OPENMP_SIMD_PRAGMA
    for (unsigned int i = 0; i < 2; ++i)
      ptr[i] = data[i];
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 16 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    _mm_stream_pd(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 2; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    for (unsigned int i = 0; i < 2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  double
  sum() const
  {
    __m128d t1 = _mm_unpackhi_pd(data, data);
    __m128d t2 = _mm_add_pd(data, t1);
    return _mm_cvtsd_f64(t2);
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
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
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_min_pd(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * Specialization for double and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const double               *in,
                              const unsigned int         *offsets,
                              VectorizedArray<double, 2> *out)
{
  const unsigned int n_chunks = n_entries / 2;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128d u0          = _mm_loadu_pd(in + 2 * i + offsets[0]);
      __m128d u1          = _mm_loadu_pd(in + 2 * i + offsets[1]);
      out[2 * i + 0].data = _mm_unpacklo_pd(u0, u1);
      out[2 * i + 1].data = _mm_unpackhi_pd(u0, u1);
    }

  // remainder loop of work that does not divide by 2
  for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 2; ++v)
      out[i][v] = in[offsets[v] + i];
}



/**
 * Specialization for double and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<double *, 2> &in,
                              VectorizedArray<double, 2>    *out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 2;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128d u0          = _mm_loadu_pd(in[0] + 2 * i);
      __m128d u1          = _mm_loadu_pd(in[1] + 2 * i);
      out[2 * i + 0].data = _mm_unpacklo_pd(u0, u1);
      out[2 * i + 1].data = _mm_unpackhi_pd(u0, u1);
    }

  for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 2; ++v)
      out[i][v] = in[v][i];
}



/**
 * Specialization for double and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 2> *in,
                               const unsigned int               *offsets,
                               double                           *out)
{
  const unsigned int n_chunks = n_entries / 2;
  if (add_into)
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out + 2 * i + offsets[0],
                        _mm_add_pd(_mm_loadu_pd(out + 2 * i + offsets[0]),
                                   res0));
          _mm_storeu_pd(out + 2 * i + offsets[1],
                        _mm_add_pd(_mm_loadu_pd(out + 2 * i + offsets[1]),
                                   res1));
        }
      // remainder loop of work that does not divide by 2
      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[offsets[v] + i] += in[i][v];
    }
  else
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out + 2 * i + offsets[0], res0);
          _mm_storeu_pd(out + 2 * i + offsets[1], res1);
        }
      // remainder loop of work that does not divide by 2
      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[offsets[v] + i] = in[i][v];
    }
}



/**
 * Specialization for double and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 2> *in,
                               std::array<double *, 2>          &out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 2;
  if (add_into)
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out[0] + 2 * i,
                        _mm_add_pd(_mm_loadu_pd(out[0] + 2 * i), res0));
          _mm_storeu_pd(out[1] + 2 * i,
                        _mm_add_pd(_mm_loadu_pd(out[1] + 2 * i), res1));
        }

      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[v][i] += in[i][v];
    }
  else
    {
      for (unsigned int i = 0; i < n_chunks; ++i)
        {
          __m128d u0   = in[2 * i + 0].data;
          __m128d u1   = in[2 * i + 1].data;
          __m128d res0 = _mm_unpacklo_pd(u0, u1);
          __m128d res1 = _mm_unpackhi_pd(u0, u1);
          _mm_storeu_pd(out[0] + 2 * i, res0);
          _mm_storeu_pd(out[1] + 2 * i, res1);
        }

      for (unsigned int i = 2 * n_chunks; i < n_entries; ++i)
        for (unsigned int v = 0; v < 2; ++v)
          out[v][i] = in[i][v];
    }
}



/**
 * Specialization for float and SSE2.
 */
template <>
class VectorizedArray<float, 4>
  : public VectorizedArrayBase<VectorizedArray<float, 4>, 4>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = float;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<float, 4>, 4>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x) &
  {
    data = _mm_set1_ps(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const float scalar) && = delete;

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  float &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm_add_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm_sub_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm_mul_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm_div_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 16 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm_loadu_ps(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 16 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm_storeu_ps(ptr, data);
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 16 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 16 == 0,
           ExcMessage("Memory not aligned"));
    _mm_stream_ps(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  float
  sum() const
  {
    __m128 t1 = _mm_movehl_ps(data, data);
    __m128 t2 = _mm_add_ps(data, t1);
    __m128 t3 = _mm_shuffle_ps(t2, t2, 1);
    __m128 t4 = _mm_add_ss(t2, t3);
    return _mm_cvtss_f32(t4);
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
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
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm_min_ps(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * Specialization for float and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int         n_entries,
                              const float               *in,
                              const unsigned int        *offsets,
                              VectorizedArray<float, 4> *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
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

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 4; ++v)
      out[i][v] = in[offsets[v] + i];
}



/**
 * Specialization for float and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int            n_entries,
                              const std::array<float *, 4> &in,
                              VectorizedArray<float, 4>    *out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m128 u0           = _mm_loadu_ps(in[0] + 4 * i);
      __m128 u1           = _mm_loadu_ps(in[1] + 4 * i);
      __m128 u2           = _mm_loadu_ps(in[2] + 4 * i);
      __m128 u3           = _mm_loadu_ps(in[3] + 4 * i);
      __m128 v0           = _mm_shuffle_ps(u0, u1, 0x44);
      __m128 v1           = _mm_shuffle_ps(u0, u1, 0xee);
      __m128 v2           = _mm_shuffle_ps(u2, u3, 0x44);
      __m128 v3           = _mm_shuffle_ps(u2, u3, 0xee);
      out[4 * i + 0].data = _mm_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm_shuffle_ps(v1, v3, 0xdd);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    for (unsigned int v = 0; v < 4; ++v)
      out[i][v] = in[v][i];
}



/**
 * Specialization for float and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 4> *in,
                               const unsigned int              *offsets,
                               float                           *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
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
      if (add_into)
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

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * Specialization for float and SSE2.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 4> *in,
                               std::array<float *, 4>          &out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
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

      if (add_into)
        {
          u0 = _mm_add_ps(_mm_loadu_ps(out[0] + 4 * i), u0);
          _mm_storeu_ps(out[0] + 4 * i, u0);
          u1 = _mm_add_ps(_mm_loadu_ps(out[1] + 4 * i), u1);
          _mm_storeu_ps(out[1] + 4 * i, u1);
          u2 = _mm_add_ps(_mm_loadu_ps(out[2] + 4 * i), u2);
          _mm_storeu_ps(out[2] + 4 * i, u2);
          u3 = _mm_add_ps(_mm_loadu_ps(out[3] + 4 * i), u3);
          _mm_storeu_ps(out[3] + 4 * i, u3);
        }
      else
        {
          _mm_storeu_ps(out[0] + 4 * i, u0);
          _mm_storeu_ps(out[1] + 4 * i, u1);
          _mm_storeu_ps(out[2] + 4 * i, u2);
          _mm_storeu_ps(out[3] + 4 * i, u3);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] = in[i][v];
}



#  endif // if DEAL_II_VECTORIZATION_WIDTH_IN_BITS > 0 && defined(__SSE2__)

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)

/**
 * Specialization of VectorizedArray class for double and AVX.
 */
template <>
class VectorizedArray<double, 4>
  : public VectorizedArrayBase<VectorizedArray<double, 4>, 4>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = double;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<double, 4>, 4>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x) &
  {
    data = _mm256_set1_pd(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const double scalar) && = delete;

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  double &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm256_add_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm256_sub_pd(data, vec.data);
#    endif
    return *this;
  }
  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm256_mul_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm256_div_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 32 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = _mm256_loadu_pd(ptr);
  }

  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm256_cvtps_pd(_mm_loadu_ps(ptr));
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 32 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    _mm256_storeu_pd(ptr, data);
  }

  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm_storeu_ps(ptr, _mm256_cvtpd_ps(data));
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 32 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 32 == 0,
           ExcMessage("Memory not aligned"));
    _mm256_stream_pd(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
#    if defined(__AVX2__) && defined(DEAL_II_USE_VECTORIZATION_GATHER)
    // unfortunately, there does not appear to be a 128 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m128 index_val =
      _mm_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m128i index = *reinterpret_cast<const __m128i *>(&index_val);

    // work around a warning with gcc-12 about an uninitialized initial state
    // for gather by starting with a zero guess, even though all lanes will be
    // overwritten
    __m256d zero = _mm256_setzero_pd();
    __m256d mask = _mm256_cmp_pd(zero, zero, _CMP_EQ_OQ);

    data = _mm256_mask_i32gather_pd(zero, base_ptr, index, mask, 8);
#    else
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
#    endif
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  double
  sum() const
  {
    VectorizedArray<double, 2> t1;
    t1.data = _mm_add_pd(this->get_lower(), this->get_upper());
    return t1.sum();
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  __m256d data;

private:
  /**
   * Extract lower half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m128d
  get_lower() const
  {
    return _mm256_castpd256_pd128(data);
  }

  /**
   * Extract upper half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m128d
  get_upper() const
  {
    return _mm256_extractf128_pd(data, 1);
  }

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
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_pd(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * Specialization for double and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const double               *in,
                              const unsigned int         *offsets,
                              VectorizedArray<double, 4> *out)
{
  const unsigned int n_chunks = n_entries / 4;
  const double      *in0      = in + offsets[0];
  const double      *in1      = in + offsets[1];
  const double      *in2      = in + offsets[2];
  const double      *in3      = in + offsets[3];

  for (unsigned int i = 0; i < n_chunks; ++i)
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

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * Specialization for double and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<double *, 4> &in,
                              VectorizedArray<double, 4>    *out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;
  const double      *in0      = in[0];
  const double      *in1      = in[1];
  const double      *in2      = in[2];
  const double      *in3      = in[3];

  for (unsigned int i = 0; i < n_chunks; ++i)
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

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * Specialization for double and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 4> *in,
                               const unsigned int               *offsets,
                               double                           *out)
{
  const unsigned int n_chunks = n_entries / 4;
  double            *out0     = out + offsets[0];
  double            *out1     = out + offsets[1];
  double            *out2     = out + offsets[2];
  double            *out3     = out + offsets[3];
  for (unsigned int i = 0; i < n_chunks; ++i)
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
      if (add_into)
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

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * Specialization for double and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 4> *in,
                               std::array<double *, 4>          &out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  double            *out0     = out[0];
  double            *out1     = out[1];
  double            *out2     = out[2];
  double            *out3     = out[3];
  for (unsigned int i = 0; i < n_chunks; ++i)
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
      if (add_into)
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

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 4; ++v)
        out[v][i] = in[i][v];
}



/**
 * Specialization for float and AVX.
 */
template <>
class VectorizedArray<float, 8>
  : public VectorizedArrayBase<VectorizedArray<float, 8>, 8>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = float;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<float, 8>, 8>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x) &
  {
    data = _mm256_set1_ps(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const float scalar) && = delete;

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  float &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m256d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm256_add_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm256_sub_ps(data, vec.data);
#    endif
    return *this;
  }
  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm256_mul_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm256_div_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 32 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm256_loadu_ps(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 32 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm256_storeu_ps(ptr, data);
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 32 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 32 == 0,
           ExcMessage("Memory not aligned"));
    _mm256_stream_ps(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
#    if defined(__AVX2__) && defined(DEAL_II_USE_VECTORIZATION_GATHER)
    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256 index_val =
      _mm256_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m256i index = *reinterpret_cast<const __m256i *>(&index_val);

    // work around a warning with gcc-12 about an uninitialized initial state
    // for gather by starting with a zero guess, even though all lanes will be
    // overwritten
    __m256 zero = _mm256_setzero_ps();
    __m256 mask = _mm256_cmp_ps(zero, zero, _CMP_EQ_OQ);

    data = _mm256_mask_i32gather_ps(zero, base_ptr, index, mask, 4);
#    else
    for (unsigned int i = 0; i < 8; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
#    endif
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    // no scatter operation in AVX/AVX2
    for (unsigned int i = 0; i < 8; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  float
  sum() const
  {
    VectorizedArray<float, 4> t1;
    t1.data = _mm_add_ps(this->get_lower(), this->get_upper());
    return t1.sum();
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  __m256 data;

private:
  /**
   * Extract lower half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m128
  get_lower() const
  {
    return _mm256_castps256_ps128(data);
  }

  /**
   * Extract upper half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m128
  get_upper() const
  {
    return _mm256_extractf128_ps(data, 1);
  }

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
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm256_min_ps(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * Specialization for float and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int         n_entries,
                              const float               *in,
                              const unsigned int        *offsets,
                              VectorizedArray<float, 8> *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      // To avoid warnings about uninitialized variables, need to initialize
      // one variable with zero before using it.
      __m256 t0, t1, t2, t3 = {};
      t0 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[0]), 0);
      t0 = _mm256_insertf128_ps(t0, _mm_loadu_ps(in + 4 * i + offsets[4]), 1);
      t1 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[1]), 0);
      t1 = _mm256_insertf128_ps(t1, _mm_loadu_ps(in + 4 * i + offsets[5]), 1);
      t2 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[2]), 0);
      t2 = _mm256_insertf128_ps(t2, _mm_loadu_ps(in + 4 * i + offsets[6]), 1);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[3]), 0);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in + 4 * i + offsets[7]), 1);

      __m256 v0           = _mm256_shuffle_ps(t0, t1, 0x44);
      __m256 v1           = _mm256_shuffle_ps(t0, t1, 0xee);
      __m256 v2           = _mm256_shuffle_ps(t2, t3, 0x44);
      __m256 v3           = _mm256_shuffle_ps(t2, t3, 0xee);
      out[4 * i + 0].data = _mm256_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm256_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm256_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm256_shuffle_ps(v1, v3, 0xdd);
    }

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * Specialization for float and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int            n_entries,
                              const std::array<float *, 8> &in,
                              VectorizedArray<float, 8>    *out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m256 t0, t1, t2, t3 = {};
      t0 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[0] + 4 * i), 0);
      t0 = _mm256_insertf128_ps(t0, _mm_loadu_ps(in[4] + 4 * i), 1);
      t1 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[1] + 4 * i), 0);
      t1 = _mm256_insertf128_ps(t1, _mm_loadu_ps(in[5] + 4 * i), 1);
      t2 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[2] + 4 * i), 0);
      t2 = _mm256_insertf128_ps(t2, _mm_loadu_ps(in[6] + 4 * i), 1);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[3] + 4 * i), 0);
      t3 = _mm256_insertf128_ps(t3, _mm_loadu_ps(in[7] + 4 * i), 1);

      __m256 v0           = _mm256_shuffle_ps(t0, t1, 0x44);
      __m256 v1           = _mm256_shuffle_ps(t0, t1, 0xee);
      __m256 v2           = _mm256_shuffle_ps(t2, t3, 0x44);
      __m256 v3           = _mm256_shuffle_ps(t2, t3, 0xee);
      out[4 * i + 0].data = _mm256_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm256_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm256_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm256_shuffle_ps(v1, v3, 0xdd);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * Specialization for float and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 8> *in,
                               const unsigned int              *offsets,
                               float                           *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
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
      if (add_into)
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

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * Specialization for float and AVX.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                       add_into,
                               const unsigned int               n_entries,
                               const VectorizedArray<float, 8> *in,
                               std::array<float *, 8>          &out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
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

      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out[0] + 4 * i), res0);
          _mm_storeu_ps(out[0] + 4 * i, res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out[1] + 4 * i), res1);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out[2] + 4 * i), res2);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out[3] + 4 * i), res3);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out[4] + 4 * i), res4);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out[5] + 4 * i), res5);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out[6] + 4 * i), res6);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out[7] + 4 * i), res7);
          _mm_storeu_ps(out[7] + 4 * i, res7);
        }
      else
        {
          _mm_storeu_ps(out[0] + 4 * i, res0);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          _mm_storeu_ps(out[7] + 4 * i, res7);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] = in[i][v];
}

#  endif

// for safety, also check that __AVX512F__ is defined in case the user manually
// set some conflicting compile flags which prevent compilation

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)

/**
 * Specialization of VectorizedArray class for double and AVX-512.
 */
template <>
class VectorizedArray<double, 8>
  : public VectorizedArrayBase<VectorizedArray<double, 8>, 8>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = double;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<double, 8>, 8>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x) &
  {
    data = _mm512_set1_pd(x);
    return *this;
  }


  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const double scalar) && = delete;

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  double &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 8);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm512_add_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm512_sub_pd(data, vec.data);
#    endif
    return *this;
  }
  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm512_mul_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm512_div_pd(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Load size() data items from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 64 bytes, as opposed
   * to casting a double address to VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = _mm512_loadu_pd(ptr);
  }

  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm512_cvtps_pd(_mm256_loadu_ps(ptr));
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 64 bytes, as opposed to casting a double address to
   * VectorizedArray<double>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    _mm512_storeu_pd(ptr, data);
  }

  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm256_storeu_ps(ptr, _mm512_cvtpd_ps(data));
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 64 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 64 == 0,
           ExcMessage("Memory not aligned"));
    _mm512_stream_pd(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
#    ifdef DEAL_II_USE_VECTORIZATION_GATHER
    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256 index_val =
      _mm256_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m256i index = *reinterpret_cast<const __m256i *>(&index_val);

    // work around a warning with gcc-12 about an uninitialized initial state
    // for gather by starting with a zero guess, even though all lanes will be
    // overwritten
    __m512d  zero = {};
    __mmask8 mask = 0xFF;

    data = _mm512_mask_i32gather_pd(zero, mask, index, base_ptr, 8);
#    else
    for (unsigned int i = 0; i < 8; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
#    endif
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
#    ifdef DEAL_II_USE_VECTORIZATION_GATHER
    for (unsigned int i = 0; i < 8; ++i)
      for (unsigned int j = i + 1; j < 8; ++j)
        Assert(offsets[i] != offsets[j],
               ExcMessage("Result of scatter undefined if two offset elements"
                          " point to the same position"));

    // unfortunately, there does not appear to be a 256 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m256 index_val =
      _mm256_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m256i index = *reinterpret_cast<const __m256i *>(&index_val);
    _mm512_i32scatter_pd(base_ptr, index, data, 8);
#    else
    for (unsigned int i = 0; i < 8; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
#    endif
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  double
  sum() const
  {
    VectorizedArray<double, 4> t1;
    t1.data = _mm256_add_pd(this->get_lower(), this->get_upper());
    return t1.sum();
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  __m512d data;

private:
  /**
   * Extract lower half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m256d
  get_lower() const
  {
    return _mm512_castpd512_pd256(data);
  }

  /**
   * Extract upper half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m256d
  get_upper() const
  {
    return _mm512_extractf64x4_pd(data, 1);
  }

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
    res.data = reinterpret_cast<__m512d>(
      _mm512_andnot_epi64(reinterpret_cast<__m512i>(mask),
                          reinterpret_cast<__m512i>(data)));
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_pd(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * Specialization for double and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const double               *in,
                              const unsigned int         *offsets,
                              VectorizedArray<double, 8> *out)
{
  // do not do full transpose because the code is long and will most
  // likely not pay off because many processors have two load units
  // (for the top 8 instructions) but only 1 permute unit (for the 8
  // shuffle/unpack instructions). rather start the transposition on the
  // vectorized array of half the size with 256 bits
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0, t1, t2, t3 = {};

      t0 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[0] + 4 * i), 0);
      t0 = _mm512_insertf64x4(t0, _mm256_loadu_pd(in + offsets[2] + 4 * i), 1);
      t1 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[1] + 4 * i), 0);
      t1 = _mm512_insertf64x4(t1, _mm256_loadu_pd(in + offsets[3] + 4 * i), 1);
      t2 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[4] + 4 * i), 0);
      t2 = _mm512_insertf64x4(t2, _mm256_loadu_pd(in + offsets[6] + 4 * i), 1);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[5] + 4 * i), 0);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in + offsets[7] + 4 * i), 1);

      __m512d v0          = _mm512_shuffle_f64x2(t0, t2, 0x88);
      __m512d v1          = _mm512_shuffle_f64x2(t0, t2, 0xdd);
      __m512d v2          = _mm512_shuffle_f64x2(t1, t3, 0x88);
      __m512d v3          = _mm512_shuffle_f64x2(t1, t3, 0xdd);
      out[4 * i + 0].data = _mm512_unpacklo_pd(v0, v2);
      out[4 * i + 1].data = _mm512_unpackhi_pd(v0, v2);
      out[4 * i + 2].data = _mm512_unpacklo_pd(v1, v3);
      out[4 * i + 3].data = _mm512_unpackhi_pd(v1, v3);
    }
  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * Specialization for double and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<double *, 8> &in,
                              VectorizedArray<double, 8>    *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0, t1, t2, t3 = {};

      t0 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[0] + 4 * i), 0);
      t0 = _mm512_insertf64x4(t0, _mm256_loadu_pd(in[2] + 4 * i), 1);
      t1 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[1] + 4 * i), 0);
      t1 = _mm512_insertf64x4(t1, _mm256_loadu_pd(in[3] + 4 * i), 1);
      t2 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[4] + 4 * i), 0);
      t2 = _mm512_insertf64x4(t2, _mm256_loadu_pd(in[6] + 4 * i), 1);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[5] + 4 * i), 0);
      t3 = _mm512_insertf64x4(t3, _mm256_loadu_pd(in[7] + 4 * i), 1);

      __m512d v0          = _mm512_shuffle_f64x2(t0, t2, 0x88);
      __m512d v1          = _mm512_shuffle_f64x2(t0, t2, 0xdd);
      __m512d v2          = _mm512_shuffle_f64x2(t1, t3, 0x88);
      __m512d v3          = _mm512_shuffle_f64x2(t1, t3, 0xdd);
      out[4 * i + 0].data = _mm512_unpacklo_pd(v0, v2);
      out[4 * i + 1].data = _mm512_unpackhi_pd(v0, v2);
      out[4 * i + 2].data = _mm512_unpacklo_pd(v1, v3);
      out[4 * i + 3].data = _mm512_unpackhi_pd(v1, v3);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * Specialization for double and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 8> *in,
                               const unsigned int               *offsets,
                               double                           *out)
{
  // as for the load, we split the store operations into 256 bit units to
  // better balance between code size, shuffle instructions, and stores
  const unsigned int n_chunks = n_entries / 4;
  __m512i mask1 = _mm512_set_epi64(0xd, 0xc, 0x5, 0x4, 0x9, 0x8, 0x1, 0x0);
  __m512i mask2 = _mm512_set_epi64(0xf, 0xe, 0x7, 0x6, 0xb, 0xa, 0x3, 0x2);
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0   = _mm512_unpacklo_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t1   = _mm512_unpackhi_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t2   = _mm512_unpacklo_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d t3   = _mm512_unpackhi_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d v0   = _mm512_permutex2var_pd(t0, mask1, t2);
      __m512d v1   = _mm512_permutex2var_pd(t0, mask2, t2);
      __m512d v2   = _mm512_permutex2var_pd(t1, mask1, t3);
      __m512d v3   = _mm512_permutex2var_pd(t1, mask2, t3);
      __m256d res0 = _mm512_extractf64x4_pd(v0, 0);
      __m256d res4 = _mm512_extractf64x4_pd(v0, 1);
      __m256d res1 = _mm512_extractf64x4_pd(v2, 0);
      __m256d res5 = _mm512_extractf64x4_pd(v2, 1);
      __m256d res2 = _mm512_extractf64x4_pd(v1, 0);
      __m256d res6 = _mm512_extractf64x4_pd(v1, 1);
      __m256d res3 = _mm512_extractf64x4_pd(v3, 0);
      __m256d res7 = _mm512_extractf64x4_pd(v3, 1);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing
      // between pointers
      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[0]), res0);
          _mm256_storeu_pd(out + 4 * i + offsets[0], res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[1]), res1);
          _mm256_storeu_pd(out + 4 * i + offsets[1], res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[2]), res2);
          _mm256_storeu_pd(out + 4 * i + offsets[2], res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[3]), res3);
          _mm256_storeu_pd(out + 4 * i + offsets[3], res3);
          res4 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[4]), res4);
          _mm256_storeu_pd(out + 4 * i + offsets[4], res4);
          res5 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[5]), res5);
          _mm256_storeu_pd(out + 4 * i + offsets[5], res5);
          res6 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[6]), res6);
          _mm256_storeu_pd(out + 4 * i + offsets[6], res6);
          res7 = _mm256_add_pd(_mm256_loadu_pd(out + 4 * i + offsets[7]), res7);
          _mm256_storeu_pd(out + 4 * i + offsets[7], res7);
        }
      else
        {
          _mm256_storeu_pd(out + 4 * i + offsets[0], res0);
          _mm256_storeu_pd(out + 4 * i + offsets[1], res1);
          _mm256_storeu_pd(out + 4 * i + offsets[2], res2);
          _mm256_storeu_pd(out + 4 * i + offsets[3], res3);
          _mm256_storeu_pd(out + 4 * i + offsets[4], res4);
          _mm256_storeu_pd(out + 4 * i + offsets[5], res5);
          _mm256_storeu_pd(out + 4 * i + offsets[6], res6);
          _mm256_storeu_pd(out + 4 * i + offsets[7], res7);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * Specialization for double and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<double, 8> *in,
                               std::array<double *, 8>          &out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  __m512i mask1 = _mm512_set_epi64(0xd, 0xc, 0x5, 0x4, 0x9, 0x8, 0x1, 0x0);
  __m512i mask2 = _mm512_set_epi64(0xf, 0xe, 0x7, 0x6, 0xb, 0xa, 0x3, 0x2);
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512d t0   = _mm512_unpacklo_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t1   = _mm512_unpackhi_pd(in[i * 4].data, in[i * 4 + 1].data);
      __m512d t2   = _mm512_unpacklo_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d t3   = _mm512_unpackhi_pd(in[i * 4 + 2].data, in[i * 4 + 3].data);
      __m512d v0   = _mm512_permutex2var_pd(t0, mask1, t2);
      __m512d v1   = _mm512_permutex2var_pd(t0, mask2, t2);
      __m512d v2   = _mm512_permutex2var_pd(t1, mask1, t3);
      __m512d v3   = _mm512_permutex2var_pd(t1, mask2, t3);
      __m256d res0 = _mm512_extractf64x4_pd(v0, 0);
      __m256d res4 = _mm512_extractf64x4_pd(v0, 1);
      __m256d res1 = _mm512_extractf64x4_pd(v2, 0);
      __m256d res5 = _mm512_extractf64x4_pd(v2, 1);
      __m256d res2 = _mm512_extractf64x4_pd(v1, 0);
      __m256d res6 = _mm512_extractf64x4_pd(v1, 1);
      __m256d res3 = _mm512_extractf64x4_pd(v3, 0);
      __m256d res7 = _mm512_extractf64x4_pd(v3, 1);

      if (add_into)
        {
          res0 = _mm256_add_pd(_mm256_loadu_pd(out[0] + 4 * i), res0);
          _mm256_storeu_pd(out[0] + 4 * i, res0);
          res1 = _mm256_add_pd(_mm256_loadu_pd(out[1] + 4 * i), res1);
          _mm256_storeu_pd(out[1] + 4 * i, res1);
          res2 = _mm256_add_pd(_mm256_loadu_pd(out[2] + 4 * i), res2);
          _mm256_storeu_pd(out[2] + 4 * i, res2);
          res3 = _mm256_add_pd(_mm256_loadu_pd(out[3] + 4 * i), res3);
          _mm256_storeu_pd(out[3] + 4 * i, res3);
          res4 = _mm256_add_pd(_mm256_loadu_pd(out[4] + 4 * i), res4);
          _mm256_storeu_pd(out[4] + 4 * i, res4);
          res5 = _mm256_add_pd(_mm256_loadu_pd(out[5] + 4 * i), res5);
          _mm256_storeu_pd(out[5] + 4 * i, res5);
          res6 = _mm256_add_pd(_mm256_loadu_pd(out[6] + 4 * i), res6);
          _mm256_storeu_pd(out[6] + 4 * i, res6);
          res7 = _mm256_add_pd(_mm256_loadu_pd(out[7] + 4 * i), res7);
          _mm256_storeu_pd(out[7] + 4 * i, res7);
        }
      else
        {
          _mm256_storeu_pd(out[0] + 4 * i, res0);
          _mm256_storeu_pd(out[1] + 4 * i, res1);
          _mm256_storeu_pd(out[2] + 4 * i, res2);
          _mm256_storeu_pd(out[3] + 4 * i, res3);
          _mm256_storeu_pd(out[4] + 4 * i, res4);
          _mm256_storeu_pd(out[5] + 4 * i, res5);
          _mm256_storeu_pd(out[6] + 4 * i, res6);
          _mm256_storeu_pd(out[7] + 4 * i, res7);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 8; ++v)
        out[v][i] = in[i][v];
}



/**
 * Specialization for float and AVX512.
 */
template <>
class VectorizedArray<float, 16>
  : public VectorizedArrayBase<VectorizedArray<float, 16>, 16>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = float;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<float, 16>, 16>(list)
  {}

  /**
   * This function can be used to set all data fields to a given scalar.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x) &
  {
    data = _mm512_set1_ps(x);
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const float scalar) && = delete;

  /**
   * Access operator.
   */
  DEAL_II_ALWAYS_INLINE
  float &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 16);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 16);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    // if the compiler supports vector arithmetic, we can simply use +=
    // operator on the given data type. this allows the compiler to combine
    // additions with multiplication (fused multiply-add) if those
    // instructions are available. Otherwise, we need to use the built-in
    // intrinsic command for __m512d
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data += vec.data;
#    else
    data = _mm512_add_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data -= vec.data;
#    else
    data = _mm512_sub_ps(data, vec.data);
#    endif
    return *this;
  }
  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data *= vec.data;
#    else
    data = _mm512_mul_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
#    ifdef DEAL_II_COMPILER_USE_VECTOR_ARITHMETICS
    data /= vec.data;
#    else
    data = _mm512_div_ps(data, vec.data);
#    endif
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address. The memory need not be aligned by 64 bytes, as opposed
   * to casting a float address to VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = _mm512_loadu_ps(ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address. The memory need not be aligned by
   * 64 bytes, as opposed to casting a float address to
   * VectorizedArray<float>*.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    _mm512_storeu_ps(ptr, data);
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   * @note Memory must be aligned by 64 bytes.
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    Assert(reinterpret_cast<std::size_t>(ptr) % 64 == 0,
           ExcMessage("Memory not aligned"));
    _mm512_stream_ps(ptr, data);
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address and with given offsets, each entry from the offset
   * providing one element of the vectorized array.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   this->operator[](v) = base_ptr[offsets[v]];
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
#    ifdef DEAL_II_USE_VECTORIZATION_GATHER
    // unfortunately, there does not appear to be a 512 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m512 index_val =
      _mm512_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m512i index = *reinterpret_cast<const __m512i *>(&index_val);

    // work around a warning with gcc-12 about an uninitialized initial state
    // for gather by starting with a zero guess, even though all lanes will be
    // overwritten
    __m512    zero = {};
    __mmask16 mask = 0xFFFF;

    data = _mm512_mask_i32gather_ps(zero, mask, index, base_ptr, 4);
#    else
    for (unsigned int i = 0; i < 16; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
#    endif
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address and the given offsets, filling the
   * elements of the vectorized array into each offset.
   *
   * This operation corresponds to the following code (but uses a more
   * efficient implementation in case the hardware allows for that):
   * @code
   * for (unsigned int v=0; v<VectorizedArray<Number>::size(); ++v)
   *   base_ptr[offsets[v]] = this->operator[](v);
   * @endcode
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
#    ifdef DEAL_II_USE_VECTORIZATION_GATHER
    for (unsigned int i = 0; i < 16; ++i)
      for (unsigned int j = i + 1; j < 16; ++j)
        Assert(offsets[i] != offsets[j],
               ExcMessage("Result of scatter undefined if two offset elements"
                          " point to the same position"));

    // unfortunately, there does not appear to be a 512 bit integer load, so
    // do it by some reinterpret casts here. this is allowed because the Intel
    // API allows aliasing between different vector types.
    const __m512 index_val =
      _mm512_loadu_ps(reinterpret_cast<const float *>(offsets));
    const __m512i index = *reinterpret_cast<const __m512i *>(&index_val);
    _mm512_i32scatter_ps(base_ptr, index, data, 4);
#    else
    for (unsigned int i = 0; i < 16; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
#    endif
  }

  /**
   * Returns sum over entries of the data field, $\sum_{i=1}^{\text{size}()}
   * this->data[i]$.
   */
  float
  sum() const
  {
    VectorizedArray<float, 8> t1;
    t1.data = _mm256_add_ps(this->get_lower(), this->get_upper());
    return t1.sum();
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  __m512 data;

private:
  /**
   * Extract lower half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m256
  get_lower() const
  {
    return _mm512_castps512_ps256(data);
  }

  /**
   * Extract upper half of data field.
   */
  DEAL_II_ALWAYS_INLINE
  __m256
  get_upper() const
  {
    return _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(data), 1));
  }

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
    res.data = reinterpret_cast<__m512>(
      _mm512_andnot_epi32(reinterpret_cast<__m512i>(mask),
                          reinterpret_cast<__m512i>(data)));
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
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
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = _mm512_min_ps(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



/**
 * Specialization for float and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int          n_entries,
                              const float                *in,
                              const unsigned int         *offsets,
                              VectorizedArray<float, 16> *out)
{
  // Similar to the double case, we perform the work on smaller entities. In
  // this case, we start from 128 bit arrays and insert them into a full 512
  // bit index. This reduces the code size and register pressure because we do
  // shuffles on 4 numbers rather than 16.
  const unsigned int n_chunks = n_entries / 4;

  // To avoid warnings about uninitialized variables, need to initialize one
  // variable to a pre-existing value in out, which will never get used in
  // the end. Keep the initialization outside the loop because of a bug in
  // gcc-9.1 which generates a "vmovapd" instruction instead of "vmovupd" in
  // case t3 is initialized to zero (inside/outside of loop), see
  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90991
  __m512 t0, t1, t2, t3;
  if (n_chunks > 0)
    t3 = out[0].data;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      t0 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[0] + 4 * i), 0);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in + offsets[4] + 4 * i), 1);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in + offsets[8] + 4 * i), 2);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in + offsets[12] + 4 * i), 3);
      t1 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[1] + 4 * i), 0);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in + offsets[5] + 4 * i), 1);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in + offsets[9] + 4 * i), 2);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in + offsets[13] + 4 * i), 3);
      t2 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[2] + 4 * i), 0);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in + offsets[6] + 4 * i), 1);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in + offsets[10] + 4 * i), 2);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in + offsets[14] + 4 * i), 3);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[3] + 4 * i), 0);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[7] + 4 * i), 1);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[11] + 4 * i), 2);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in + offsets[15] + 4 * i), 3);

      __m512 v0 = _mm512_shuffle_ps(t0, t1, 0x44);
      __m512 v1 = _mm512_shuffle_ps(t0, t1, 0xee);
      __m512 v2 = _mm512_shuffle_ps(t2, t3, 0x44);
      __m512 v3 = _mm512_shuffle_ps(t2, t3, 0xee);

      out[4 * i + 0].data = _mm512_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm512_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm512_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm512_shuffle_ps(v1, v3, 0xdd);
    }

  // remainder loop of work that does not divide by 4
  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    out[i].gather(in + i, offsets);
}



/**
 * Specialization for float and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_load_and_transpose(const unsigned int             n_entries,
                              const std::array<float *, 16> &in,
                              VectorizedArray<float, 16>    *out)
{
  // see the comments in the vectorized_load_and_transpose above

  const unsigned int n_chunks = n_entries / 4;

  __m512 t0, t1, t2, t3;
  if (n_chunks > 0)
    t3 = out[0].data;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      t0 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[0] + 4 * i), 0);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in[4] + 4 * i), 1);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in[8] + 4 * i), 2);
      t0 = _mm512_insertf32x4(t0, _mm_loadu_ps(in[12] + 4 * i), 3);
      t1 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[1] + 4 * i), 0);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in[5] + 4 * i), 1);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in[9] + 4 * i), 2);
      t1 = _mm512_insertf32x4(t1, _mm_loadu_ps(in[13] + 4 * i), 3);
      t2 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[2] + 4 * i), 0);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in[6] + 4 * i), 1);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in[10] + 4 * i), 2);
      t2 = _mm512_insertf32x4(t2, _mm_loadu_ps(in[14] + 4 * i), 3);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[3] + 4 * i), 0);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[7] + 4 * i), 1);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[11] + 4 * i), 2);
      t3 = _mm512_insertf32x4(t3, _mm_loadu_ps(in[15] + 4 * i), 3);

      __m512 v0 = _mm512_shuffle_ps(t0, t1, 0x44);
      __m512 v1 = _mm512_shuffle_ps(t0, t1, 0xee);
      __m512 v2 = _mm512_shuffle_ps(t2, t3, 0x44);
      __m512 v3 = _mm512_shuffle_ps(t2, t3, 0xee);

      out[4 * i + 0].data = _mm512_shuffle_ps(v0, v2, 0x88);
      out[4 * i + 1].data = _mm512_shuffle_ps(v0, v2, 0xdd);
      out[4 * i + 2].data = _mm512_shuffle_ps(v1, v3, 0x88);
      out[4 * i + 3].data = _mm512_shuffle_ps(v1, v3, 0xdd);
    }

  for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
    gather(out[i], in, i);
}



/**
 * Specialization for float and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<float, 16> *in,
                               const unsigned int               *offsets,
                               float                            *out)
{
  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512 t0 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0x44);
      __m512 t1 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0xee);
      __m512 t2 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0x44);
      __m512 t3 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0xee);
      __m512 u0 = _mm512_shuffle_ps(t0, t2, 0x88);
      __m512 u1 = _mm512_shuffle_ps(t0, t2, 0xdd);
      __m512 u2 = _mm512_shuffle_ps(t1, t3, 0x88);
      __m512 u3 = _mm512_shuffle_ps(t1, t3, 0xdd);

      __m128 res0  = _mm512_extractf32x4_ps(u0, 0);
      __m128 res4  = _mm512_extractf32x4_ps(u0, 1);
      __m128 res8  = _mm512_extractf32x4_ps(u0, 2);
      __m128 res12 = _mm512_extractf32x4_ps(u0, 3);
      __m128 res1  = _mm512_extractf32x4_ps(u1, 0);
      __m128 res5  = _mm512_extractf32x4_ps(u1, 1);
      __m128 res9  = _mm512_extractf32x4_ps(u1, 2);
      __m128 res13 = _mm512_extractf32x4_ps(u1, 3);
      __m128 res2  = _mm512_extractf32x4_ps(u2, 0);
      __m128 res6  = _mm512_extractf32x4_ps(u2, 1);
      __m128 res10 = _mm512_extractf32x4_ps(u2, 2);
      __m128 res14 = _mm512_extractf32x4_ps(u2, 3);
      __m128 res3  = _mm512_extractf32x4_ps(u3, 0);
      __m128 res7  = _mm512_extractf32x4_ps(u3, 1);
      __m128 res11 = _mm512_extractf32x4_ps(u3, 2);
      __m128 res15 = _mm512_extractf32x4_ps(u3, 3);

      // Cannot use the same store instructions in both paths of the 'if'
      // because the compiler cannot know that there is no aliasing between
      // pointers
      if (add_into)
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
          res8 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[8]), res8);
          _mm_storeu_ps(out + 4 * i + offsets[8], res8);
          res9 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[9]), res9);
          _mm_storeu_ps(out + 4 * i + offsets[9], res9);
          res10 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[10]), res10);
          _mm_storeu_ps(out + 4 * i + offsets[10], res10);
          res11 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[11]), res11);
          _mm_storeu_ps(out + 4 * i + offsets[11], res11);
          res12 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[12]), res12);
          _mm_storeu_ps(out + 4 * i + offsets[12], res12);
          res13 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[13]), res13);
          _mm_storeu_ps(out + 4 * i + offsets[13], res13);
          res14 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[14]), res14);
          _mm_storeu_ps(out + 4 * i + offsets[14], res14);
          res15 = _mm_add_ps(_mm_loadu_ps(out + 4 * i + offsets[15]), res15);
          _mm_storeu_ps(out + 4 * i + offsets[15], res15);
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
          _mm_storeu_ps(out + 4 * i + offsets[8], res8);
          _mm_storeu_ps(out + 4 * i + offsets[9], res9);
          _mm_storeu_ps(out + 4 * i + offsets[10], res10);
          _mm_storeu_ps(out + 4 * i + offsets[11], res11);
          _mm_storeu_ps(out + 4 * i + offsets[12], res12);
          _mm_storeu_ps(out + 4 * i + offsets[13], res13);
          _mm_storeu_ps(out + 4 * i + offsets[14], res14);
          _mm_storeu_ps(out + 4 * i + offsets[15], res15);
        }
    }

  // remainder loop of work that does not divide by 4
  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[offsets[v] + i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[offsets[v] + i] = in[i][v];
}



/**
 * Specialization for float and AVX-512.
 */
template <>
inline DEAL_II_ALWAYS_INLINE void
vectorized_transpose_and_store(const bool                        add_into,
                               const unsigned int                n_entries,
                               const VectorizedArray<float, 16> *in,
                               std::array<float *, 16>          &out)
{
  // see the comments in the vectorized_transpose_and_store above

  const unsigned int n_chunks = n_entries / 4;
  for (unsigned int i = 0; i < n_chunks; ++i)
    {
      __m512 t0 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0x44);
      __m512 t1 = _mm512_shuffle_ps(in[4 * i].data, in[1 + 4 * i].data, 0xee);
      __m512 t2 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0x44);
      __m512 t3 =
        _mm512_shuffle_ps(in[2 + 4 * i].data, in[3 + 4 * i].data, 0xee);
      __m512 u0 = _mm512_shuffle_ps(t0, t2, 0x88);
      __m512 u1 = _mm512_shuffle_ps(t0, t2, 0xdd);
      __m512 u2 = _mm512_shuffle_ps(t1, t3, 0x88);
      __m512 u3 = _mm512_shuffle_ps(t1, t3, 0xdd);

      __m128 res0  = _mm512_extractf32x4_ps(u0, 0);
      __m128 res4  = _mm512_extractf32x4_ps(u0, 1);
      __m128 res8  = _mm512_extractf32x4_ps(u0, 2);
      __m128 res12 = _mm512_extractf32x4_ps(u0, 3);
      __m128 res1  = _mm512_extractf32x4_ps(u1, 0);
      __m128 res5  = _mm512_extractf32x4_ps(u1, 1);
      __m128 res9  = _mm512_extractf32x4_ps(u1, 2);
      __m128 res13 = _mm512_extractf32x4_ps(u1, 3);
      __m128 res2  = _mm512_extractf32x4_ps(u2, 0);
      __m128 res6  = _mm512_extractf32x4_ps(u2, 1);
      __m128 res10 = _mm512_extractf32x4_ps(u2, 2);
      __m128 res14 = _mm512_extractf32x4_ps(u2, 3);
      __m128 res3  = _mm512_extractf32x4_ps(u3, 0);
      __m128 res7  = _mm512_extractf32x4_ps(u3, 1);
      __m128 res11 = _mm512_extractf32x4_ps(u3, 2);
      __m128 res15 = _mm512_extractf32x4_ps(u3, 3);

      if (add_into)
        {
          res0 = _mm_add_ps(_mm_loadu_ps(out[0] + 4 * i), res0);
          _mm_storeu_ps(out[0] + 4 * i, res0);
          res1 = _mm_add_ps(_mm_loadu_ps(out[1] + 4 * i), res1);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          res2 = _mm_add_ps(_mm_loadu_ps(out[2] + 4 * i), res2);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          res3 = _mm_add_ps(_mm_loadu_ps(out[3] + 4 * i), res3);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          res4 = _mm_add_ps(_mm_loadu_ps(out[4] + 4 * i), res4);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          res5 = _mm_add_ps(_mm_loadu_ps(out[5] + 4 * i), res5);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          res6 = _mm_add_ps(_mm_loadu_ps(out[6] + 4 * i), res6);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          res7 = _mm_add_ps(_mm_loadu_ps(out[7] + 4 * i), res7);
          _mm_storeu_ps(out[7] + 4 * i, res7);
          res8 = _mm_add_ps(_mm_loadu_ps(out[8] + 4 * i), res8);
          _mm_storeu_ps(out[8] + 4 * i, res8);
          res9 = _mm_add_ps(_mm_loadu_ps(out[9] + 4 * i), res9);
          _mm_storeu_ps(out[9] + 4 * i, res9);
          res10 = _mm_add_ps(_mm_loadu_ps(out[10] + 4 * i), res10);
          _mm_storeu_ps(out[10] + 4 * i, res10);
          res11 = _mm_add_ps(_mm_loadu_ps(out[11] + 4 * i), res11);
          _mm_storeu_ps(out[11] + 4 * i, res11);
          res12 = _mm_add_ps(_mm_loadu_ps(out[12] + 4 * i), res12);
          _mm_storeu_ps(out[12] + 4 * i, res12);
          res13 = _mm_add_ps(_mm_loadu_ps(out[13] + 4 * i), res13);
          _mm_storeu_ps(out[13] + 4 * i, res13);
          res14 = _mm_add_ps(_mm_loadu_ps(out[14] + 4 * i), res14);
          _mm_storeu_ps(out[14] + 4 * i, res14);
          res15 = _mm_add_ps(_mm_loadu_ps(out[15] + 4 * i), res15);
          _mm_storeu_ps(out[15] + 4 * i, res15);
        }
      else
        {
          _mm_storeu_ps(out[0] + 4 * i, res0);
          _mm_storeu_ps(out[1] + 4 * i, res1);
          _mm_storeu_ps(out[2] + 4 * i, res2);
          _mm_storeu_ps(out[3] + 4 * i, res3);
          _mm_storeu_ps(out[4] + 4 * i, res4);
          _mm_storeu_ps(out[5] + 4 * i, res5);
          _mm_storeu_ps(out[6] + 4 * i, res6);
          _mm_storeu_ps(out[7] + 4 * i, res7);
          _mm_storeu_ps(out[8] + 4 * i, res8);
          _mm_storeu_ps(out[9] + 4 * i, res9);
          _mm_storeu_ps(out[10] + 4 * i, res10);
          _mm_storeu_ps(out[11] + 4 * i, res11);
          _mm_storeu_ps(out[12] + 4 * i, res12);
          _mm_storeu_ps(out[13] + 4 * i, res13);
          _mm_storeu_ps(out[14] + 4 * i, res14);
          _mm_storeu_ps(out[15] + 4 * i, res15);
        }
    }

  if (add_into)
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[v][i] += in[i][v];
  else
    for (unsigned int i = 4 * n_chunks; i < n_entries; ++i)
      for (unsigned int v = 0; v < 16; ++v)
        out[v][i] = in[i][v];
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ALTIVEC__) && \
    defined(__VSX__)

template <>
class VectorizedArray<double, 2>
  : public VectorizedArrayBase<VectorizedArray<double, 2>, 2>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = double;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const double scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<double, 2>, 2>(list)
  {}

  /**
   * This function assigns a scalar to the current object.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const double x) &
  {
    data = vec_splats(x);

    // Some compilers believe that vec_splats sets 'x', but that's not true.
    // They then warn about setting a variable and not using it. Suppress the
    // warning by "using" the variable:
    (void)x;
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const double scalar) && = delete;

  /**
   * Access operator. The component must be either 0 or 1.
   */
  DEAL_II_ALWAYS_INLINE
  double &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<double *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const double &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 2);
    return *(reinterpret_cast<const double *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data = vec_add(data, vec.data);
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data = vec_sub(data, vec.data);
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data = vec_mul(data, vec.data);
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data = vec_div(data, vec.data);
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const double *ptr)
  {
    data = vec_vsx_ld(0, ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(double *ptr) const
  {
    vec_vsx_st(data, 0, ptr);
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(double *ptr) const
  {
    store(ptr);
  }

  /**
   * @copydoc VectorizedArray<Number>::gather()
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const double *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 2; ++i)
      *(reinterpret_cast<double *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * @copydoc VectorizedArray<Number>::scatter
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, double *base_ptr) const
  {
    for (unsigned int i = 0; i < 2; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const double *>(&data) + i);
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  __vector double data;

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
    res.data = vec_sqrt(data);
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
    res.data = vec_abs(data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_max(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_min(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};



template <>
class VectorizedArray<float, 4>
  : public VectorizedArrayBase<VectorizedArray<float, 4>, 4>
{
public:
  /**
   * This gives the type of the array elements.
   */
  using value_type = float;

  /**
   * Record the fact that the given specialization of VectorizedArray is
   * indeed implemented.
   */
  static constexpr bool is_implemented = true;

  /**
   * Default empty constructor, leaving the data in an uninitialized state
   * similar to float/double.
   */
  VectorizedArray() = default;

  /**
   * Construct an array with the given scalar broadcast to all lanes.
   */
  VectorizedArray(const float scalar)
  {
    this->operator=(scalar);
  }

  /**
   * Construct an array with the given initializer list.
   */
  template <typename U>
  VectorizedArray(const std::initializer_list<U> &list)
    : VectorizedArrayBase<VectorizedArray<float, 4>, 4>(list)
  {}

  /**
   * This function assigns a scalar to the current object.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator=(const float x) &
  {
    data = vec_splats(x);

    // Some compilers believe that vec_splats sets 'x', but that's not true.
    // They then warn about setting a variable and not using it. Suppress the
    // warning by "using" the variable:
    (void)x;
    return *this;
  }

  /**
   * Assign a scalar to the current object. This overload is used for
   * rvalue references; because it does not make sense to assign
   * something to a temporary, the function is deleted.
   */
  VectorizedArray &
  operator=(const float scalar) && = delete;

  /**
   * Access operator. The component must be between 0 and 3.
   */
  DEAL_II_ALWAYS_INLINE
  float &
  operator[](const unsigned int comp)
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<float *>(&data) + comp);
  }

  /**
   * Constant access operator.
   */
  DEAL_II_ALWAYS_INLINE
  const float &
  operator[](const unsigned int comp) const
  {
    AssertIndexRange(comp, 4);
    return *(reinterpret_cast<const float *>(&data) + comp);
  }

  /**
   * Element-wise addition of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator+=(const VectorizedArray &vec)
  {
    data = vec_add(data, vec.data);
    return *this;
  }

  /**
   * Element-wise subtraction of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator-=(const VectorizedArray &vec)
  {
    data = vec_sub(data, vec.data);
    return *this;
  }

  /**
   * Element-wise multiplication of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator*=(const VectorizedArray &vec)
  {
    data = vec_mul(data, vec.data);
    return *this;
  }

  /**
   * Element-wise division of two arrays of numbers.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray &
  operator/=(const VectorizedArray &vec)
  {
    data = vec_div(data, vec.data);
    return *this;
  }

  /**
   * Load @p size() from memory into the calling class, starting at
   * the given address.
   */
  DEAL_II_ALWAYS_INLINE
  void
  load(const float *ptr)
  {
    data = vec_vsx_ld(0, ptr);
  }

  /**
   * Write the content of the calling class into memory in form of @p
   * size() to the given address.
   */
  DEAL_II_ALWAYS_INLINE
  void
  store(float *ptr) const
  {
    vec_vsx_st(data, 0, ptr);
  }

  /**
   * @copydoc VectorizedArray<Number>::streaming_store()
   */
  DEAL_II_ALWAYS_INLINE
  void
  streaming_store(float *ptr) const
  {
    store(ptr);
  }

  /**
   * @copydoc VectorizedArray<Number>::gather()
   */
  DEAL_II_ALWAYS_INLINE
  void
  gather(const float *base_ptr, const unsigned int *offsets)
  {
    for (unsigned int i = 0; i < 4; ++i)
      *(reinterpret_cast<float *>(&data) + i) = base_ptr[offsets[i]];
  }

  /**
   * @copydoc VectorizedArray<Number>::scatter
   */
  DEAL_II_ALWAYS_INLINE
  void
  scatter(const unsigned int *offsets, float *base_ptr) const
  {
    for (unsigned int i = 0; i < 4; ++i)
      base_ptr[offsets[i]] = *(reinterpret_cast<const float *>(&data) + i);
  }

  /**
   * Actual data field. To be consistent with the standard layout type and to
   * enable interaction with external SIMD functionality, this member is
   * declared public.
   */
  __vector float data;

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
    res.data = vec_sqrt(data);
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
    res.data = vec_abs(data);
    return res;
  }

  /**
   * Return the component-wise maximum of this field and another one. Not for
   * use in user code. Use max(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_max(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_max(data, other.data);
    return res;
  }

  /**
   * Return the component-wise minimum of this field and another one. Not for
   * use in user code. Use min(x,y) instead.
   */
  DEAL_II_ALWAYS_INLINE
  VectorizedArray
  get_min(const VectorizedArray &other) const
  {
    VectorizedArray res;
    res.data = vec_min(data, other.data);
    return res;
  }

  // Make a few functions friends.
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::sqrt(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::abs(const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::max(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
  template <typename Number2, std::size_t width2>
  friend VectorizedArray<Number2, width2>
  std::min(const VectorizedArray<Number2, width2> &,
           const VectorizedArray<Number2, width2> &);
};

#  endif // if DEAL_II_VECTORIZATION_LEVEL >=1 && defined(__ALTIVEC__) &&
         // defined(__VSX__)


#endif // DOXYGEN



/**
 * @name Arithmetic operations with VectorizedArray
 * @{
 */

/**
 * Relational operator == for VectorizedArray
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE bool
operator==(const VectorizedArray<Number, width> &lhs,
           const VectorizedArray<Number, width> &rhs)
{
  for (unsigned int i = 0; i < VectorizedArray<Number, width>::size(); ++i)
    if (lhs[i] != rhs[i])
      return false;

  return true;
}


/**
 * Addition of two vectorized arrays with operator +.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator+(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp += v;
}

/**
 * Subtraction of two vectorized arrays with operator -.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator-(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp -= v;
}

/**
 * Multiplication of two vectorized arrays with operator *.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator*(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp *= v;
}

/**
 * Division of two vectorized arrays with operator /.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator/(const VectorizedArray<Number, width> &u,
          const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp /= v;
}

/**
 * Addition of a scalar (expanded to a vectorized array with @p
 * size() equal entries) and a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator+(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp += v;
}

/**
 * Addition of a scalar (expanded to a vectorized array with @p
 * size() equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator+(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = u;
  return tmp += v;
}

/**
 * Addition of a vectorized array and a scalar (expanded to a vectorized array
 * with @p size() equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator+(const VectorizedArray<Number, width> &v, const Number &u)
{
  return u + v;
}

/**
 * Addition of a vectorized array and a scalar (expanded to a vectorized array
 * with @p size() equal entries) in case the scalar is a double
 * (needed in order to be able to write simple code with constants that are
 * usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator+(const VectorizedArray<float, width> &v, const double u)
{
  return u + v;
}

/**
 * Subtraction of a vectorized array from a scalar (expanded to a vectorized
 * array with @p size() equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator-(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp -= v;
}

/**
 * Subtraction of a vectorized array from a scalar (expanded to a vectorized
 * array with @p size() equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator-(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return tmp -= v;
}

/**
 * Subtraction of a scalar (expanded to a vectorized array with @p
 * size() equal entries) from a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator-(const VectorizedArray<Number, width> &v, const Number &u)
{
  VectorizedArray<Number, width> tmp = u;
  return v - tmp;
}

/**
 * Subtraction of a scalar (expanded to a vectorized array with @p
 * size() equal entries) from a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator-(const VectorizedArray<float, width> &v, const double u)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return v - tmp;
}

/**
 * Multiplication of a scalar (expanded to a vectorized array with @p
 * size() equal entries) and a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator*(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp *= v;
}

/**
 * Multiplication of a scalar (expanded to a vectorized array with @p
 * size() equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator*(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return tmp *= v;
}

/**
 * Multiplication of a vectorized array and a scalar (expanded to a vectorized
 * array with @p size() equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator*(const VectorizedArray<Number, width> &v, const Number &u)
{
  return u * v;
}

/**
 * Multiplication of a vectorized array and a scalar (expanded to a vectorized
 * array with @p size() equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator*(const VectorizedArray<float, width> &v, const double u)
{
  return u * v;
}

/**
 * Quotient between a scalar (expanded to a vectorized array with @p
 * size() equal entries) and a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator/(const Number &u, const VectorizedArray<Number, width> &v)
{
  VectorizedArray<Number, width> tmp = u;
  return tmp /= v;
}

/**
 * Quotient between a scalar (expanded to a vectorized array with @p
 * size() equal entries) and a vectorized array in case the scalar
 * is a double (needed in order to be able to write simple code with constants
 * that are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator/(const double u, const VectorizedArray<float, width> &v)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return tmp /= v;
}

/**
 * Quotient between a vectorized array and a scalar (expanded to a vectorized
 * array with @p size() equal entries).
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
operator/(const VectorizedArray<Number, width> &v, const Number &u)
{
  VectorizedArray<Number, width> tmp = u;
  return v / tmp;
}

/**
 * Quotient between a vectorized array and a scalar (expanded to a vectorized
 * array with @p size() equal entries) in case the scalar is a
 * double (needed in order to be able to write simple code with constants that
 * are usually double numbers).
 *
 * @relatesalso VectorizedArray
 */
template <std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<float, width>
operator/(const VectorizedArray<float, width> &v, const double u)
{
  VectorizedArray<float, width> tmp = static_cast<float>(u);
  return v / tmp;
}

/**
 * Unary operator + on a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator+(const VectorizedArray<Number, width> &u)
{
  return u;
}

/**
 * Unary operator - on a vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline DEAL_II_ALWAYS_INLINE VectorizedArray<Number, width>
                             operator-(const VectorizedArray<Number, width> &u)
{
  // to get a negative sign, subtract the input from zero (could also
  // multiply by -1, but this one is slightly simpler)
  return VectorizedArray<Number, width>() - u;
}

/**
 * Output operator for vectorized array.
 *
 * @relatesalso VectorizedArray
 */
template <typename Number, std::size_t width>
inline std::ostream &
operator<<(std::ostream &out, const VectorizedArray<Number, width> &p)
{
  constexpr unsigned int n = VectorizedArray<Number, width>::size();
  for (unsigned int i = 0; i < n - 1; ++i)
    out << p[i] << ' ';
  out << p[n - 1];

  return out;
}

/** @} */

/**
 * @name Ternary operations on VectorizedArray
 * @{
 */

/**
 * enum class encoding binary operations for a component-wise comparison of
 * VectorizedArray data types.
 *
 * @note In case of SIMD vecorization (sse, avx, av512) we select the
 * corresponding ordered, non-signalling (<code>OQ</code>) variants.
 */
enum class SIMDComparison : int
{
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)
  equal                 = _CMP_EQ_OQ,
  not_equal             = _CMP_NEQ_OQ,
  less_than             = _CMP_LT_OQ,
  less_than_or_equal    = _CMP_LE_OQ,
  greater_than          = _CMP_GT_OQ,
  greater_than_or_equal = _CMP_GE_OQ
#else
  equal,
  not_equal,
  less_than,
  less_than_or_equal,
  greater_than,
  greater_than_or_equal
#endif
};


/**
 * Computes the vectorized equivalent of the following ternary operation:
 * @code
 *   (left OP right) ? true_value : false_value
 * @endcode
 * where <code>OP</code> is a binary operator (such as <code>=</code>,
 * <code>!=</code>, <code><</code>, <code><=</code>, <code>></code>, and
 * <code>>=</code>).
 *
 * Such a computational idiom is useful as an alternative to branching
 * whenever the control flow itself would depend on (computed) data. For
 * example, in case of a scalar data type the statement
 * <code>(left < right) ? true_value : false_value</code>
 * could have been also implemented using an <code>if</code>-statement:
 * @code
 * if (left < right)
 *     result = true_value;
 * else
 *     result = false_value;
 * @endcode
 * This, however, is fundamentally impossible in case of vectorization
 * because different decisions will be necessary on different vector entries
 * (lanes) and
 * the first variant (based on a ternary operator) has to be used instead:
 * @code
 *   result = compare_and_apply_mask<SIMDComparison::less_than>
 *     (left, right, true_value, false_value);
 * @endcode
 * Some more illustrative examples (that are less efficient than the
 * dedicated <code>std::max</code> and <code>std::abs</code> overloads):
 * @code
 *   VectorizedArray<double> left;
 *   VectorizedArray<double> right;
 *
 *   // std::max
 *   const auto maximum = compare_and_apply_mask<SIMDComparison::greater_than>
 *     (left, right, left, right);
 *
 *   // std::abs
 *   const auto absolute = compare_and_apply_mask<SIMDComparison::less_than>
 *     (left, VectorizedArray<double>(0.), -left, left);
 * @endcode
 *
 * More precisely, this function first computes a (boolean) mask that is
 * the result of a binary operator <code>OP</code> applied to all elements
 * of the VectorizedArray arguments @p left and @p right. The mask is then
 * used to either select the corresponding component of @p true_value (if
 * the binary operation equates to true), or @p false_value. The binary
 * operator is encoded via the SIMDComparison template argument
 * @p predicate.
 *
 * In order to ease with generic programming approaches, the function
 * provides overloads for all VectorizedArray<Number> variants as well as
 * generic POD types such as double and float.
 *
 * @note For this function to work the binary operation has to be encoded
 * via a SIMDComparison template argument @p predicate. Depending on it
 * appropriate low-level machine instructions are generated replacing the
 * call to compare_and_apply_mask. This also explains why @p predicate is a
 * compile-time constant template parameter and not a constant function
 * argument. In order to be able to emit the correct low-level instruction,
 * the compiler has to know the comparison at compile time.
 */
template <SIMDComparison predicate, typename Number>
DEAL_II_ALWAYS_INLINE inline Number
compare_and_apply_mask(const Number &left,
                       const Number &right,
                       const Number &true_value,
                       const Number &false_value)
{
  bool mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = (left == right);
        break;
      case SIMDComparison::not_equal:
        mask = (left != right);
        break;
      case SIMDComparison::less_than:
        mask = (left < right);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = (left <= right);
        break;
      case SIMDComparison::greater_than:
        mask = (left > right);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = (left >= right);
        break;
    }

  return mask ? true_value : false_value;
}


/**
 * Specialization of above function for the non-vectorized
 * VectorizedArray<Number, 1> variant.
 */
template <SIMDComparison predicate, typename Number>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<Number, 1>
compare_and_apply_mask(const VectorizedArray<Number, 1> &left,
                       const VectorizedArray<Number, 1> &right,
                       const VectorizedArray<Number, 1> &true_value,
                       const VectorizedArray<Number, 1> &false_value)
{
  VectorizedArray<Number, 1> result;
  result.data = compare_and_apply_mask<predicate, Number>(left.data,
                                                          right.data,
                                                          true_value.data,
                                                          false_value.data);
  return result;
}

/** @} */

#ifndef DOXYGEN
#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 16>
compare_and_apply_mask(const VectorizedArray<float, 16> &left,
                       const VectorizedArray<float, 16> &right,
                       const VectorizedArray<float, 16> &true_values,
                       const VectorizedArray<float, 16> &false_values)
{
  const __mmask16 mask =
    _mm512_cmp_ps_mask(left.data, right.data, static_cast<int>(predicate));
  VectorizedArray<float, 16> result;
  result.data = _mm512_mask_mov_ps(false_values.data, mask, true_values.data);
  return result;
}



template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 8>
compare_and_apply_mask(const VectorizedArray<double, 8> &left,
                       const VectorizedArray<double, 8> &right,
                       const VectorizedArray<double, 8> &true_values,
                       const VectorizedArray<double, 8> &false_values)
{
  const __mmask16 mask =
    _mm512_cmp_pd_mask(left.data, right.data, static_cast<int>(predicate));
  VectorizedArray<double, 8> result;
  result.data = _mm512_mask_mov_pd(false_values.data, mask, true_values.data);
  return result;
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 8>
compare_and_apply_mask(const VectorizedArray<float, 8> &left,
                       const VectorizedArray<float, 8> &right,
                       const VectorizedArray<float, 8> &true_values,
                       const VectorizedArray<float, 8> &false_values)
{
  const auto mask =
    _mm256_cmp_ps(left.data, right.data, static_cast<int>(predicate));

  VectorizedArray<float, 8> result;
  result.data = _mm256_blendv_ps(false_values.data, true_values.data, mask);
  return result;
}


template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 4>
compare_and_apply_mask(const VectorizedArray<double, 4> &left,
                       const VectorizedArray<double, 4> &right,
                       const VectorizedArray<double, 4> &true_values,
                       const VectorizedArray<double, 4> &false_values)
{
  const auto mask =
    _mm256_cmp_pd(left.data, right.data, static_cast<int>(predicate));

  VectorizedArray<double, 4> result;
  result.data = _mm256_blendv_pd(false_values.data, true_values.data, mask);
  return result;
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 4>
compare_and_apply_mask(const VectorizedArray<float, 4> &left,
                       const VectorizedArray<float, 4> &right,
                       const VectorizedArray<float, 4> &true_values,
                       const VectorizedArray<float, 4> &false_values)
{
  __m128 mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = _mm_cmpeq_ps(left.data, right.data);
        break;
      case SIMDComparison::not_equal:
        mask = _mm_cmpneq_ps(left.data, right.data);
        break;
      case SIMDComparison::less_than:
        mask = _mm_cmplt_ps(left.data, right.data);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = _mm_cmple_ps(left.data, right.data);
        break;
      case SIMDComparison::greater_than:
        mask = _mm_cmpgt_ps(left.data, right.data);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = _mm_cmpge_ps(left.data, right.data);
        break;
    }

  VectorizedArray<float, 4> result;
  result.data = _mm_or_ps(_mm_and_ps(mask, true_values.data),
                          _mm_andnot_ps(mask, false_values.data));

  return result;
}


template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 2>
compare_and_apply_mask(const VectorizedArray<double, 2> &left,
                       const VectorizedArray<double, 2> &right,
                       const VectorizedArray<double, 2> &true_values,
                       const VectorizedArray<double, 2> &false_values)
{
  __m128d mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = _mm_cmpeq_pd(left.data, right.data);
        break;
      case SIMDComparison::not_equal:
        mask = _mm_cmpneq_pd(left.data, right.data);
        break;
      case SIMDComparison::less_than:
        mask = _mm_cmplt_pd(left.data, right.data);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = _mm_cmple_pd(left.data, right.data);
        break;
      case SIMDComparison::greater_than:
        mask = _mm_cmpgt_pd(left.data, right.data);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = _mm_cmpge_pd(left.data, right.data);
        break;
    }

  VectorizedArray<double, 2> result;
  result.data = _mm_or_pd(_mm_and_pd(mask, true_values.data),
                          _mm_andnot_pd(mask, false_values.data));

  return result;
}

#  endif

#  if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ARM_NEON)

template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<float, 4>
compare_and_apply_mask(const VectorizedArray<float, 4> &left,
                       const VectorizedArray<float, 4> &right,
                       const VectorizedArray<float, 4> &true_values,
                       const VectorizedArray<float, 4> &false_values)
{
  uint32x4_t mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = vceqq_f32(left.data, right.data);
        break;
      case SIMDComparison::not_equal:
        mask = vmvnq_u32(vceqq_f32(left.data, right.data));
        break;
      case SIMDComparison::less_than:
        mask = vcltq_f32(left.data, right.data);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = vcleq_f32(left.data, right.data);
        break;
      case SIMDComparison::greater_than:
        mask = vcgtq_f32(left.data, right.data);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = vcgeq_f32(left.data, right.data);
        break;
    }

  VectorizedArray<float, 4> result;
  result.data = vreinterpretq_f32_u32(vorrq_u32(
    vandq_u32(mask, vreinterpretq_u32_f32(true_values.data)),
    vandq_u32(vmvnq_u32(mask), vreinterpretq_u32_f32(false_values.data))));

  return result;
}


template <SIMDComparison predicate>
DEAL_II_ALWAYS_INLINE inline VectorizedArray<double, 2>
compare_and_apply_mask(const VectorizedArray<double, 2> &left,
                       const VectorizedArray<double, 2> &right,
                       const VectorizedArray<double, 2> &true_values,
                       const VectorizedArray<double, 2> &false_values)
{
  uint64x2_t mask;
  switch (predicate)
    {
      case SIMDComparison::equal:
        mask = vceqq_f64(left.data, right.data);
        break;
      case SIMDComparison::not_equal:
        mask = vreinterpretq_u64_u32(
          vmvnq_u32(vreinterpretq_u32_u64(vceqq_f64(left.data, right.data))));
        break;
      case SIMDComparison::less_than:
        mask = vcltq_f64(left.data, right.data);
        break;
      case SIMDComparison::less_than_or_equal:
        mask = vcleq_f64(left.data, right.data);
        break;
      case SIMDComparison::greater_than:
        mask = vcgtq_f64(left.data, right.data);
        break;
      case SIMDComparison::greater_than_or_equal:
        mask = vcgeq_f64(left.data, right.data);
        break;
    }

  VectorizedArray<double, 2> result;
  result.data = vreinterpretq_f64_u64(vorrq_u64(
    vandq_u64(mask, vreinterpretq_u64_f64(true_values.data)),
    vandq_u64(vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(mask))),
              vreinterpretq_u64_f64(false_values.data))));

  return result;
}

#  endif
#endif // DOXYGEN


namespace internal
{
  template <typename T>
  struct VectorizedArrayTrait
  {
    /**
     * Define scalar value type.
     */
    using value_type = T;

    /**
     * Return the width of template type.
     */
    static constexpr std::size_t
    width()
    {
      return 1;
    }

    /**
     * Define vectorized value type for internal vectorization.
     */
    using vectorized_value_type = VectorizedArray<T>;

    /**
     * Return a stride which defines how often the template type T fits into
     * the vectorized_value_type. This is useful to write vectorized templated
     * code where the internal computation is vectorized and the user
     * interface is optionally scalar or also vectorized.
     */
    static constexpr std::size_t
    stride()
    {
      return vectorized_value_type::size();
    }

    /**
     * Get a reference to scalar value (on lane 0).
     */
    static value_type &
    get(value_type &value, unsigned int c)
    {
      AssertIndexRange(c, 1);
      (void)c;

      return value;
    }

    /**
     * Get a read-only reference to scalar value (on lane 0).
     */
    static const value_type &
    get(const value_type &value, unsigned int c)
    {
      AssertIndexRange(c, 1);
      (void)c;

      return value;
    }

    /**
     * Get a reference to scalar value on lane c from a vectorized values field.
     */
    static value_type &
    get_from_vectorized(vectorized_value_type &values, unsigned int c)
    {
      AssertIndexRange(c, stride());

      return values[c];
    }

    /**
     * Get a read-only reference to scalar value on lane c from a vectorized
     * values field.
     */
    static const value_type &
    get_from_vectorized(const vectorized_value_type &values, unsigned int c)
    {
      AssertIndexRange(c, stride());

      return values[c];
    }
  };

  template <typename T, std::size_t width_>
  struct VectorizedArrayTrait<VectorizedArray<T, width_>>
  {
    /**
     * Define scalar value type.
     */
    using value_type = T;

    /**
     * Return the width of template type.
     */
    static constexpr std::size_t
    width()
    {
      return width_;
    }

    /**
     * Define vectorized value type for internal vectorization.
     */
    using vectorized_value_type = VectorizedArray<T, width_>;

    /**
     * Return a stride which defines how often the template type
     * VectorizedArray<T, width_> fits into the vectorized value type. This is
     * useful to write vectorized templated code where the internal computation
     * is vectorized and the user interface is optionally scalar or also
     * vectorized.
     */
    static constexpr std::size_t
    stride()
    {
      return 1;
    }

    /**
     * Get a reference to scalar value on lane c.
     */
    static value_type &
    get(vectorized_value_type &values, unsigned int c)
    {
      AssertIndexRange(c, width_);

      return values[c];
    }

    /**
     * Get a read-only reference to scalar value on lane c.
     */
    static const value_type &
    get(const vectorized_value_type &values, unsigned int c)
    {
      AssertIndexRange(c, width_);

      return values[c];
    }

    /**
     * Get a reference to vectorized values from a vectorized values field.
     */
    static vectorized_value_type &
    get_from_vectorized(vectorized_value_type &values, unsigned int c)
    {
      (void)c;
      AssertIndexRange(c, stride());

      return values;
    }

    /**
     * Get a read-only reference to vectorized values from a vectorized values
     * field.
     */
    static const vectorized_value_type &
    get_from_vectorized(const vectorized_value_type &values, unsigned int c)
    {
      (void)c;
      AssertIndexRange(c, stride());

      return values;
    }
  };
} // namespace internal


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
   * sin(x[VectorizedArray::size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  sin(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::sin(x[i]);
    return out;
  }



  /**
   * Compute the cosine of a vectorized data field. The result is returned as
   * vectorized array in the form <tt>{cos(x[0]), cos(x[1]), ...,
   * cos(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  cos(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::cos(x[i]);
    return out;
  }



  /**
   * Compute the tangent of a vectorized data field. The result is returned
   * as vectorized array in the form <tt>{tan(x[0]), tan(x[1]), ...,
   * tan(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  tan(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::tan(x[i]);
    return out;
  }



  /**
   * Compute the arc cosine of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{acos(x[0]), acos(x[1]), ...,
   * acos(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  acos(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::acos(x[i]);
    return out;
  }



  /**
   * Compute the arc sine of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{asin(x[0]), asin(x[1]), ...,
   * asin(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  asin(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::asin(x[i]);
    return out;
  }



  /**
   * Compute the arc tangent of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{atan(x[0]), atan(x[1]), ...,
   * atan(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  atan(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::atan(x[i]);
    return out;
  }



  /**
   * Compute the hyperbolic cosine of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{cosh(x[0]), cosh(x[1]), ...,
   * cosh(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  cosh(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::cosh(x[i]);
    return out;
  }



  /**
   * Compute the hyperbolic sine of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{sinh(x[0]), sinh(x[1]), ...,
   * sinh(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  sinh(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::sinh(x[i]);
    return out;
  }



  /**
   * Compute the hyperbolic tangent of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{tanh(x[0]), tanh(x[1]), ...,
   * tanh(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  tanh(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::tanh(x[i]);
    return out;
  }



  /**
   * Compute the area hyperbolic cosine of a vectorized data field. The result
   * is returned as vectorized array in the form <tt>{acosh(x[0]), acosh(x[1]),
   * ..., acosh(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  acosh(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::acosh(x[i]);
    return out;
  }



  /**
   * Compute the area hyperbolic sine of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{asinh(x[0]), asinh(x[1]),
   * ..., asinh(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  asinh(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::asinh(x[i]);
    return out;
  }



  /**
   * Compute the area hyperbolic tangent of a vectorized data field. The result
   * is returned as vectorized array in the form <tt>{atanh(x[0]), atanh(x[1]),
   * ..., atanh(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  atanh(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::atanh(x[i]);
    return out;
  }



  /**
   * Compute the exponential of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{exp(x[0]), exp(x[1]), ...,
   * exp(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  exp(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::exp(x[i]);
    return out;
  }



  /**
   * Compute the natural logarithm of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{log(x[0]), log(x[1]), ...,
   * log(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  log(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::log(x[i]);
    return out;
  }



  /**
   * Compute the square root of a vectorized data field. The result is
   * returned as vectorized array in the form <tt>{sqrt(x[0]), sqrt(x[1]),
   * ..., sqrt(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  sqrt(const ::dealii::VectorizedArray<Number, width> &x)
  {
    return x.get_sqrt();
  }



  /**
   * Raises the given number @p x to the power @p p for a vectorized data
   * field. The result is returned as vectorized array in the form
   * <tt>{pow(x[0],p), pow(x[1],p), ..., pow(x[size()-1],p)}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  pow(const ::dealii::VectorizedArray<Number, width> &x, const Number p)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::pow(x[i], p);
    return out;
  }



  /**
   * Raises the given number @p x to the power @p p for a vectorized data
   * field. The result is returned as vectorized array in the form
   * <tt>{pow(x[0],p[0]), pow(x[1],p[1]), ...,
   * pow(x[size()-1],p[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  pow(const ::dealii::VectorizedArray<Number, width> &x,
      const ::dealii::VectorizedArray<Number, width> &p)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
         ++i)
      out[i] = std::pow(x[i], p[i]);
    return out;
  }



  /**
   * Compute the absolute value (modulus) of a vectorized data field. The
   * result is returned as vectorized array in the form <tt>{abs(x[0]),
   * abs(x[1]), ..., abs(x[size()-1])}</tt>.
   *
   * @relatesalso VectorizedArray
   */
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  abs(const ::dealii::VectorizedArray<Number, width> &x)
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
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  max(const ::dealii::VectorizedArray<Number, width> &x,
      const ::dealii::VectorizedArray<Number, width> &y)
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
  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  min(const ::dealii::VectorizedArray<Number, width> &x,
      const ::dealii::VectorizedArray<Number, width> &y)
  {
    return x.get_min(y);
  }



  /**
   * Iterator traits for VectorizedArrayIterator.
   */
  template <class T>
  struct iterator_traits<dealii::VectorizedArrayIterator<T>>
  {
#ifdef DEAL_II_HAVE_CXX20
    using iterator_category = contiguous_iterator_tag;
#else
    using iterator_category = random_access_iterator_tag;
#endif
    using value_type      = typename T::value_type;
    using difference_type = std::ptrdiff_t;
  };

} // namespace std

#endif
