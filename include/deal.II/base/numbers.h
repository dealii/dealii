// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2019 by the deal.II authors
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

#ifndef dealii_numbers_h
#define dealii_numbers_h


#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  include <cuComplex.h>
#endif

#include <cmath>
#include <complex>
#include <cstddef>
#include <type_traits>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  define DEAL_II_CUDA_HOST_DEV __host__ __device__
#else
#  define DEAL_II_CUDA_HOST_DEV
#endif

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * A helper class specifying the maximal vector length of VectorizedArray
   * for a specified data type Number for the given processor architecture and
   * optimization level.
   *
   * The value of the maximal vector length is used as default template
   * argument in VectorizedArray, such that VectorizedArray<Number> is
   * equivalent to VectorizedArray<Number,
   * VectorizedArrayWidthSpecifier<Number>::max_width>.
   *
   * @note This class is the default implementation for data types for which
   * no vectorization is supported.
   *
   * @tparam Number The underlying data type for which one wants to find out
   *   the maximal length of hardware supported vectors.
   *
   * @author Peter Munch, 2019
   */
  template <typename Number>
  struct VectorizedArrayWidthSpecifier
  {
    /**
     * Maximal vector length of VectorizedArray for an arbitrary type.
     */
    constexpr static unsigned int max_width = 1;
  };

  /**
   * A helper class specifying the maximal vector length of VectorizedArray
   * for the data type `double` for the given processor architecture and
   * optimization level. For a detailed description of supported maximal vector
   * lengths, see the the documentation of VectorizedArray.
   *
   * @author Peter Munch, 2019
   */
  template <>
  struct VectorizedArrayWidthSpecifier<double>
  {
    /**
     * Maximal vector length of VectorizedArray for double.
     */
    constexpr static unsigned int max_width =
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
      8;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
      4;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
      2;
#else
      1;
#endif
  };

  /**
   * A helper class specifying the maximal vector length of VectorizedArray
   * for the data type `float` for the given processor architecture and
   * optimization level. For a detailed description of supported maximal vector
   * lengths, see the the documentation of VectorizedArray.
   *
   * @author Peter Munch, 2019
   */
  template <>
  struct VectorizedArrayWidthSpecifier<float>
  {
    /**
     * Maximal vector length of VectorizedArray for float.
     */
    constexpr static unsigned int max_width =
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__ALTIVEC__)
      4;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512 && defined(__AVX512F__)
      16;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256 && defined(__AVX__)
      8;
#elif DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128 && defined(__SSE2__)
      4;
#else
      1;
#endif
  };


} // namespace internal

// forward declarations to support abs or sqrt operations on VectorizedArray
#ifndef DOXYGEN
template <typename Number,
          std::size_t width =
            internal::VectorizedArrayWidthSpecifier<Number>::max_width>
class VectorizedArray;
template <typename T>
struct EnableIfScalar;
#endif

DEAL_II_NAMESPACE_CLOSE

// Declare / Import auto-differentiable math functions in(to) standard
// namespace before numbers::NumberTraits is defined
#ifdef DEAL_II_WITH_ADOLC
#  include <deal.II/differentiation/ad/adolc_math.h>

#  include <adolc/adouble.h> // Taped double
#endif
// Ideally we'd like to #include <deal.II/differentiation/ad/sacado_math.h>
// but header indirectly references numbers.h. We therefore simply
// import the whole Sacado header at this point to get the math
// functions imported into the standard namespace.
#ifdef DEAL_II_TRILINOS_WITH_SACADO
#  include <Sacado.hpp>
#endif

namespace std
{
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  sqrt(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  abs(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  max(const ::dealii::VectorizedArray<Number, width> &,
      const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, std::size_t width>
  DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number, width>
  min(const ::dealii::VectorizedArray<Number, width> &,
      const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  pow(const ::dealii::VectorizedArray<Number, width> &, const Number p);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  sin(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  cos(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  tan(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  exp(const ::dealii::VectorizedArray<Number, width> &);
  template <typename Number, size_t width>
  ::dealii::VectorizedArray<Number, width>
  log(const ::dealii::VectorizedArray<Number, width> &);
} // namespace std

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for the declaration of universal constants. Since the
 * availability in <tt>math.h</tt> is not always guaranteed, we put them here.
 * Since this file is included by <tt>base/config.h</tt>, they are available
 * to the whole library.
 *
 * The constants defined here are a subset of the <tt>M_XXX</tt> constants
 * sometimes declared in the system include file <tt>math.h</tt>, but without
 * the prefix <tt>M_</tt>.
 *
 * In addition to that, we declare  <tt>invalid_unsigned_int</tt> to be the
 * largest unsigned integer representable; this value is widely used in the
 * library as a marker for an invalid index, an invalid size of an array, and
 * similar purposes.
 */
namespace numbers
{
  /**
   * e
   */
  static constexpr double E = 2.7182818284590452354;

  /**
   * log_2 e
   */
  static constexpr double LOG2E = 1.4426950408889634074;

  /**
   * log_10 e
   */
  static constexpr double LOG10E = 0.43429448190325182765;

  /**
   * log_e 2
   */
  static constexpr double LN2 = 0.69314718055994530942;

  /**
   * log_e 10
   */
  static constexpr double LN10 = 2.30258509299404568402;

  /**
   * pi
   */
  static constexpr double PI = 3.14159265358979323846;

  /**
   * pi/2
   */
  static constexpr double PI_2 = 1.57079632679489661923;

  /**
   * pi/4
   */
  static constexpr double PI_4 = 0.78539816339744830962;

  /**
   * sqrt(2)
   */
  static constexpr double SQRT2 = 1.41421356237309504880;

  /**
   * 1/sqrt(2)
   */
  static constexpr double SQRT1_2 = 0.70710678118654752440;

  /**
   * Check whether the given type can be used in CUDA device code.
   * If not, DEAL_II_CUDA_HOST_DEV needs to be disabled for functions
   * that use this type.
   */
  template <typename Number, typename = void>
  struct is_cuda_compatible : std::true_type
  {};

  /**
   * std::complex cannot be used in CUDA device code.
   */
  template <typename Number>
  struct is_cuda_compatible<std::complex<Number>, void> : std::false_type
  {};

  /**
   * Check whether a value is not a number.
   *
   * This function uses either <code>std::isnan</code>, <code>isnan</code>, or
   * <code>_isnan</code>, whichever is available on the system and returns the
   * result.
   *
   * If none of the functions detecting NaN is available, this function
   * returns false.
   *
   * @deprecated This function has been deprecated in favor of the C++11
   * function <code>std::isnan</code>.
   */
  DEAL_II_DEPRECATED
  bool
  is_nan(const double x);

  /**
   * Return @p true if the given value is a finite floating point number, i.e.
   * is neither plus or minus infinity nor NaN (not a number).
   *
   * Note that the argument type of this function is <code>double</code>. In
   * other words, if you give a very large number of type <code>long
   * double</code>, this function may return <code>false</code> even if the
   * number is finite with respect to type <code>long double</code>.
   */
  bool
  is_finite(const double x);

  /**
   * Return @p true if real and imaginary parts of the given complex number
   * are finite.
   */
  bool
  is_finite(const std::complex<double> &x);

  /**
   * Return @p true if real and imaginary parts of the given complex number
   * are finite.
   */
  bool
  is_finite(const std::complex<float> &x);

  /**
   * Return @p true if real and imaginary parts of the given complex number
   * are finite.
   *
   * Again may not work correctly if real or imaginary parts are very large
   * numbers that are infinite in terms of <code>double</code>, but finite
   * with respect to <code>long double</code>.
   */
  bool
  is_finite(const std::complex<long double> &x);

  /**
   * Return whether two numbers are equal to one another.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note This function expects that @p value_2 is castable to the type
   * of @p value_1.
   */
  template <typename Number1, typename Number2>
  constexpr bool
  values_are_equal(const Number1 &value_1, const Number2 &value_2);

  /**
   * Return whether two numbers are not equal to one another.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note This function expects that @p value_2 is castable to the type
   * of @p value_1.
   */
  template <typename Number1, typename Number2>
  bool
  values_are_not_equal(const Number1 &value_1, const Number2 &value_2);

  /**
   * Return whether or not a value is equal to zero.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   */
  template <typename Number>
  constexpr bool
  value_is_zero(const Number &value);

  /**
   * Return whether @p value_1 is less than that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note This function expects that @p value_2 is castable to the type
   * of @p value_1.
   */
  template <typename Number1, typename Number2>
  bool
  value_is_less_than(const Number1 &value_1, const Number2 &value_2);

  /**
   * Return whether @p value_1 is less than or equal to that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note This function expects that @p value_2 is castable to the type
   * of @p value_1.
   */
  template <typename Number1, typename Number2>
  bool
  value_is_less_than_or_equal_to(const Number1 &value_1,
                                 const Number2 &value_2);



  /**
   * Return whether @p value_1 is greater than that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note This function expects that @p value_2 is castable to the type
   * of @p value_1.
   */
  template <typename Number1, typename Number2>
  bool
  value_is_greater_than(const Number1 &value_1, const Number2 &value_2);

  /**
   * Return whether @p value_1 is greater than or equal to that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note This function expects that @p value_2 is castable to the type
   * of @p value_1.
   */
  template <typename Number1, typename Number2>
  bool
  value_is_greater_than_or_equal_to(const Number1 &value_1,
                                    const Number2 &value_2);

  /**
   * A structure that, together with its partial specializations
   * NumberTraits<std::complex<number> >, provides traits and member functions
   * that make it possible to write templates that work on both real number
   * types and complex number types. This template is mostly used to implement
   * linear algebra classes such as vectors and matrices that work for both
   * real and complex numbers.
   *
   * @author Wolfgang Bangerth, 2007
   */
  template <typename number>
  struct NumberTraits
  {
    /**
     * A flag that specifies whether the template type given to this class is
     * complex or real. Since the general template is selected for non-complex
     * types, the answer is <code>false</code>.
     */
    static constexpr bool is_complex = false;

    /**
     * For this data type, alias the corresponding real type. Since the
     * general template is selected for all data types that are not
     * specializations of std::complex<T>, the underlying type must be real-
     * values, so the real_type is equal to the underlying type.
     */
    using real_type = number;

    /**
     * Return the complex-conjugate of the given number. Since the general
     * template is selected if number is not a complex data type, this
     * function simply returns the given number.
     *
     * @note This function can also be used in CUDA device code.
     */
    static constexpr DEAL_II_CUDA_HOST_DEV const number &
                                                 conjugate(const number &x);

    /**
     * Return the square of the absolute value of the given number. Since the
     * general template is chosen for types not equal to std::complex, this
     * function simply returns the square of the given number.
     *
     * @note If the template type can be used in CUDA device code, the same holds true
     * for this function.
     */
    template <typename Dummy = number>
    static constexpr DEAL_II_CUDA_HOST_DEV
      typename std::enable_if<std::is_same<Dummy, number>::value &&
                                is_cuda_compatible<Dummy>::value,
                              real_type>::type
      abs_square(const number &x);

    template <typename Dummy = number>
    static constexpr
      typename std::enable_if<std::is_same<Dummy, number>::value &&
                                !is_cuda_compatible<Dummy>::value,
                              real_type>::type
      abs_square(const number &x);

    /**
     * Return the absolute value of a number.
     */
    static real_type
    abs(const number &x);
  };


  /**
   * Specialization of the general NumberTraits class that provides the
   * relevant information if the underlying data type is std::complex<T>.
   *
   * @author Wolfgang Bangerth, 2007
   */
  template <typename number>
  struct NumberTraits<std::complex<number>>
  {
    /**
     * A flag that specifies whether the template type given to this class is
     * complex or real. Since this specialization of the general template is
     * selected for complex types, the answer is <code>true</code>.
     */
    static constexpr bool is_complex = true;

    /**
     * For this data type, alias the corresponding real type. Since this
     * specialization of the template is selected for number types
     * std::complex<T>, the real type is equal to the type used to store the
     * two components of the complex number.
     */
    using real_type = number;

    /**
     * Return the complex-conjugate of the given number.
     */
    static constexpr std::complex<number>
    conjugate(const std::complex<number> &x);

    /**
     * Return the square of the absolute value of the given number. Since this
     * specialization of the general template is chosen for types equal to
     * std::complex, this function returns the product of a number and its
     * complex conjugate.
     */
    static constexpr real_type
    abs_square(const std::complex<number> &x);


    /**
     * Return the absolute value of a complex number.
     */
    static real_type
    abs(const std::complex<number> &x);
  };

  // --------------- inline and template functions ---------------- //

  inline bool
  is_nan(const double x)
  {
    return std::isnan(x);
  }



  inline bool
  is_finite(const double x)
  {
    return std::isfinite(x);
  }



  inline bool
  is_finite(const std::complex<double> &x)
  {
    // Check complex numbers for infinity
    // by testing real and imaginary part
    return (is_finite(x.real()) && is_finite(x.imag()));
  }



  inline bool
  is_finite(const std::complex<float> &x)
  {
    // Check complex numbers for infinity
    // by testing real and imaginary part
    return (is_finite(x.real()) && is_finite(x.imag()));
  }



  inline bool
  is_finite(const std::complex<long double> &x)
  {
    // Same for std::complex<long double>
    return (is_finite(x.real()) && is_finite(x.imag()));
  }


  template <typename number>
  constexpr DEAL_II_CUDA_HOST_DEV const number &
                                        NumberTraits<number>::conjugate(const number &x)
  {
    return x;
  }



  template <typename number>
  template <typename Dummy>
  constexpr DEAL_II_CUDA_HOST_DEV
    typename std::enable_if<std::is_same<Dummy, number>::value &&
                              is_cuda_compatible<Dummy>::value,
                            typename NumberTraits<number>::real_type>::type
    NumberTraits<number>::abs_square(const number &x)
  {
    return x * x;
  }



  template <typename number>
  template <typename Dummy>
  constexpr
    typename std::enable_if<std::is_same<Dummy, number>::value &&
                              !is_cuda_compatible<Dummy>::value,
                            typename NumberTraits<number>::real_type>::type
    NumberTraits<number>::abs_square(const number &x)
  {
    return x * x;
  }



  template <typename number>
  typename NumberTraits<number>::real_type
  NumberTraits<number>::abs(const number &x)
  {
    return std::abs(x);
  }



  template <typename number>
  constexpr std::complex<number>
  NumberTraits<std::complex<number>>::conjugate(const std::complex<number> &x)
  {
    return std::conj(x);
  }



  template <typename number>
  typename NumberTraits<std::complex<number>>::real_type
  NumberTraits<std::complex<number>>::abs(const std::complex<number> &x)
  {
    return std::abs(x);
  }



  template <typename number>
  constexpr typename NumberTraits<std::complex<number>>::real_type
  NumberTraits<std::complex<number>>::abs_square(const std::complex<number> &x)
  {
    return std::norm(x);
  }

} // namespace numbers


// Forward declarations
namespace Differentiation
{
  namespace AD
  {
    namespace internal
    {
      // Defined in differentiation/ad/ad_number_traits.h
      template <typename T>
      struct NumberType;
    } // namespace internal

    // Defined in differentiation/ad/ad_number_traits.h
    template <typename NumberType>
    struct is_ad_number;
  } // namespace AD
} // namespace Differentiation


namespace internal
{
  /**
   * A test to see if it is possible to convert one number
   * type to the other.
   */
  template <typename From, typename To>
  struct is_explicitly_convertible
  {
    // Source: https://stackoverflow.com/a/16944130
  private:
    template <typename T>
    static void f(T);

    template <typename F, typename T>
    static constexpr auto
    test(int) -> decltype(f(static_cast<T>(std::declval<F>())), true)
    {
      return true;
    }

    template <typename F, typename T>
    static constexpr auto
    test(...) -> bool
    {
      return false;
    }

  public:
    static bool const value = test<From, To>(0);
  };

  /*
   * The structs below are needed to convert between some special number types.
   * Also see tensor.h for another specialization.
   */
  template <typename T>
  struct NumberType
  {
    static constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV const T &
                                                                       value(const T &t)
    {
      return t;
    }

    // Below are generic functions that allows an overload for any
    // type U that is transformable to type T. This is particularly
    // useful when needing to cast exotic number types
    // (e.g. auto-differentiable or symbolic numbers) to a floating
    // point one, such as might  happen when converting between tensor
    // types.

    // Type T is constructible from F.
    template <typename F>
    static constexpr DEAL_II_ALWAYS_INLINE DEAL_II_CUDA_HOST_DEV T
                                                                 value(const F &f,
                                                                       typename std::enable_if<
            !std::is_same<typename std::decay<T>::type,
                          typename std::decay<F>::type>::value &&
            std::is_constructible<T, F>::value>::type * = nullptr)
    {
      return T(f);
    }

    // Type T is explicitly convertible (but not constructible) from F.
    template <typename F>
    static constexpr DEAL_II_ALWAYS_INLINE T
                                           value(const F &f,
                                                 typename std::enable_if<
            !std::is_same<typename std::decay<T>::type,
                          typename std::decay<F>::type>::value &&
            !std::is_constructible<T, F>::value &&
            is_explicitly_convertible<const F, T>::value>::type * = nullptr)
    {
      return static_cast<T>(f);
    }

    // Sacado doesn't provide any conversion operators, so we have
    // to extract the value and perform further conversions from there.
    // To be safe, we extend this to other possible AD numbers that
    // might fall into the same category.
    template <typename F>
    static T
    value(const F &f,
          typename std::enable_if<
            !std::is_same<typename std::decay<T>::type,
                          typename std::decay<F>::type>::value &&
            !std::is_constructible<T, F>::value &&
            !is_explicitly_convertible<const F, T>::value &&
            Differentiation::AD::is_ad_number<F>::value>::type * = nullptr)
    {
      return Differentiation::AD::internal::NumberType<T>::value(f);
    }
  };

  template <typename T>
  struct NumberType<std::complex<T>>
  {
    static constexpr const std::complex<T> &
    value(const std::complex<T> &t)
    {
      return t;
    }

    static constexpr std::complex<T>
    value(const T &t)
    {
      return std::complex<T>(t);
    }

    // Facilitate cast from complex<double> to complex<float>
    template <typename U>
    static constexpr std::complex<T>
    value(const std::complex<U> &t)
    {
      return std::complex<T>(NumberType<T>::value(t.real()),
                             NumberType<T>::value(t.imag()));
    }
  };

#ifdef DEAL_II_COMPILER_CUDA_AWARE
  template <>
  struct NumberType<cuComplex>
  {
    static cuComplex
    value(const float t)
    {
      return make_cuComplex(t, 0.f);
    }
  };

  template <>
  struct NumberType<cuDoubleComplex>
  {
    static cuDoubleComplex
    value(const double t)
    {
      return make_cuDoubleComplex(t, 0.);
    }
  };
#endif
} // namespace internal

namespace numbers
{
#ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

  /**
   * Return whether two numbers are equal to one another. For intricate data
   * types (e.g. some automatically differentiable numbers), this function
   * returns only whether the scalar values stored by the input values are
   * equal.
   *
   * @note When ADOL-C is compiled with the "advanced branching" feature, then
   * this specialization is only intended for use in assertions and
   * other code paths that do not affect the end result of a computation.
   */
  // Defined in differentiation/ad/adolc_number_types.cc
  bool
  values_are_equal(const adouble &value_1, const adouble &value_2);


  /**
   * Return whether two numbers are equal to one another. For intricate data
   * types (e.g. some automatically differentiable numbers), this function
   * returns only whether the scalar values stored by the input values are
   * equal.
   *
   * @note When ADOL-C is compiled with the "advanced branching" feature, then
   * this specialization is only intended for use in assertions and
   * other code paths that do not affect the end result of a computation.
   */
  template <typename Number>
  bool
  values_are_equal(const adouble &value_1, const Number &value_2)
  {
    // Use the specialized definition for two ADOL-C taped types
    return values_are_equal(value_1,
                            internal::NumberType<adouble>::value(value_2));
  }


  /**
   * Return whether two numbers are equal to one another. For intricate data
   * types (e.g. some automatically differentiable numbers), this function
   * returns only whether the scalar values stored by the input values are
   * equal.
   *
   * @note When ADOL-C is compiled with the "advanced branching" feature, then
   * this specialization is only intended for use in assertions and
   * other code paths that do not affect the end result of a computation.
   */
  template <typename Number>
  bool
  values_are_equal(const Number &value_1, const adouble &value_2)
  {
    // Use the above definition
    return values_are_equal(value_2, value_1);
  }

  /**
   * Return whether @p value_1 is less than that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note When ADOL-C is compiled with the "advanced branching" feature, then
   * this specialization is only intended for use in assertions and
   * other code paths that do not affect the end result of a computation.
   */
  // Defined in differentiation/ad/adolc_number_types.cc
  bool
  value_is_less_than(const adouble &value_1, const adouble &value_2);


  /**
   * Return whether @p value_1 is less than that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note When ADOL-C is compiled with the "advanced branching" feature, then
   * this specialization is only intended for use in assertions and
   * other code paths that do not affect the end result of a computation.
   */
  template <typename Number>
  bool
  value_is_less_than(const adouble &value_1, const Number &value_2)
  {
    // Use the specialized definition for two ADOL-C taped types
    return value_is_less_than(value_1,
                              internal::NumberType<adouble>::value(value_2));
  }


  /**
   * Return whether @p value_1 is less than that of @p value_2.
   *
   * For intricate data types (e.g. some automatically differentiable numbers),
   * this function returns the result of the comparison of scalar values stored
   * by the input arguments.
   *
   * @note When ADOL-C is compiled with the "advanced branching" feature, then
   * this specialization is only intended for use in assertions and
   * other code paths that do not affect the end result of a computation.
   */
  template <typename Number>
  bool
  value_is_less_than(const Number &value_1, const adouble &value_2)
  {
    // Use the specialized definition for two ADOL-C taped types
    return value_is_less_than(internal::NumberType<adouble>::value(value_1),
                              value_2);
  }

#endif


  template <typename Number1, typename Number2>
  constexpr bool
  values_are_equal(const Number1 &value_1, const Number2 &value_2)
  {
    return (value_1 == internal::NumberType<Number1>::value(value_2));
  }


  template <typename Number1, typename Number2>
  inline bool
  values_are_not_equal(const Number1 &value_1, const Number2 &value_2)
  {
    return !(values_are_equal(value_1, value_2));
  }


  template <typename Number>
  constexpr bool
  value_is_zero(const Number &value)
  {
    return values_are_equal(value, 0.0);
  }


  template <typename Number1, typename Number2>
  inline bool
  value_is_less_than(const Number1 &value_1, const Number2 &value_2)
  {
    return (value_1 < internal::NumberType<Number1>::value(value_2));
  }


  template <typename Number1, typename Number2>
  inline bool
  value_is_less_than_or_equal_to(const Number1 &value_1, const Number2 &value_2)
  {
    return (value_is_less_than(value_1, value_2) ||
            values_are_equal(value_1, value_2));
  }


  template <typename Number1, typename Number2>
  bool
  value_is_greater_than(const Number1 &value_1, const Number2 &value_2)
  {
    return !(value_is_less_than_or_equal_to(value_1, value_2));
  }


  template <typename Number1, typename Number2>
  inline bool
  value_is_greater_than_or_equal_to(const Number1 &value_1,
                                    const Number2 &value_2)
  {
    return !(value_is_less_than(value_1, value_2));
  }
} // namespace numbers

DEAL_II_NAMESPACE_CLOSE

#endif
