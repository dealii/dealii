// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2017 by the deal.II authors
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

#ifndef dealii_numbers_h
#define dealii_numbers_h


#include <deal.II/base/config.h>
#include <deal.II/base/types.h>

#include <cmath>
#include <cstdlib>
#include <complex>

#ifdef DEAL_II_WITH_CUDA
#  include <cuda_runtime_api.h>
#  define DEAL_II_CUDA_HOST_DEV __host__ __device__
#else
#  define DEAL_II_CUDA_HOST_DEV
#endif

DEAL_II_NAMESPACE_OPEN

// forward declarations to support abs or sqrt operations on VectorizedArray
template <typename Number> class VectorizedArray;
template <typename T> struct EnableIfScalar;

DEAL_II_NAMESPACE_CLOSE

namespace std
{
  template <typename Number> DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number>
  sqrt(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number>
  abs(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number>
  max(const ::dealii::VectorizedArray<Number> &, const ::dealii::VectorizedArray<Number> &);
  template <typename Number> DEAL_II_ALWAYS_INLINE ::dealii::VectorizedArray<Number>
  min (const ::dealii::VectorizedArray<Number> &, const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  pow(const ::dealii::VectorizedArray<Number> &, const Number p);
  template <typename Number> ::dealii::VectorizedArray<Number>
  sin(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  cos(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  tan(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  exp(const ::dealii::VectorizedArray<Number> &);
  template <typename Number> ::dealii::VectorizedArray<Number>
  log(const ::dealii::VectorizedArray<Number> &);
}

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
  static const double  E       = 2.7182818284590452354;

  /**
   * log_2 e
   */
  static const double  LOG2E   = 1.4426950408889634074;

  /**
   * log_10 e
   */
  static const double  LOG10E  = 0.43429448190325182765;

  /**
   * log_e 2
   */
  static const double  LN2     = 0.69314718055994530942;

  /**
   * log_e 10
   */
  static const double  LN10    = 2.30258509299404568402;

  /**
   * pi
   */
  static const double  PI      = 3.14159265358979323846;

  /**
   * pi/2
   */
  static const double  PI_2    = 1.57079632679489661923;

  /**
   * pi/4
   */
  static const double  PI_4    = 0.78539816339744830962;

  /**
   * sqrt(2)
   */
  static const double  SQRT2   = 1.41421356237309504880;

  /**
   * 1/sqrt(2)
   */
  static const double  SQRT1_2 = 0.70710678118654752440;

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
  bool is_nan (const double x);

  /**
   * Return @p true if the given value is a finite floating point number, i.e.
   * is neither plus or minus infinity nor NaN (not a number).
   *
   * Note that the argument type of this function is <code>double</code>. In
   * other words, if you give a very large number of type <code>long
   * double</code>, this function may return <code>false</code> even if the
   * number is finite with respect to type <code>long double</code>.
   */
  bool is_finite (const double x);

  /**
   * Return @p true if real and imaginary parts of the given complex number
   * are finite.
   */
  bool is_finite (const std::complex<double> &x);

  /**
   * Return @p true if real and imaginary parts of the given complex number
   * are finite.
   */
  bool is_finite (const std::complex<float> &x);

  /**
   * Return @p true if real and imaginary parts of the given complex number
   * are finite.
   *
   * Again may not work correctly if real or imaginary parts are very large
   * numbers that are infinite in terms of <code>double</code>, but finite
   * with respect to <code>long double</code>.
   */
  bool is_finite (const std::complex<long double> &x);

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
    static const bool is_complex = false;

    /**
     * For this data type, typedef the corresponding real type. Since the
     * general template is selected for all data types that are not
     * specializations of std::complex<T>, the underlying type must be real-
     * values, so the real_type is equal to the underlying type.
     */
    typedef number real_type;

    /**
     * Return the complex-conjugate of the given number. Since the general
     * template is selected if number is not a complex data type, this
     * function simply returns the given number.
     */
    static
    DEAL_II_CUDA_HOST_DEV
    const number &conjugate (const number &x);

    /**
     * Return the square of the absolute value of the given number. Since the
     * general template is chosen for types not equal to std::complex, this
     * function simply returns the square of the given number.
     *
     * @ingroup CUDAWrappers
     */
    static
    DEAL_II_CUDA_HOST_DEV
    real_type abs_square (const number &x);

    /**
     * Return the absolute value of a number.
     */
    static
    real_type abs (const number &x);
  };


  /**
   * Specialization of the general NumberTraits class that provides the
   * relevant information if the underlying data type is std::complex<T>.
   *
   * @author Wolfgang Bangerth, 2007
   */
  template <typename number>
  struct NumberTraits<std::complex<number> >
  {
    /**
     * A flag that specifies whether the template type given to this class is
     * complex or real. Since this specialization of the general template is
     * selected for complex types, the answer is <code>true</code>.
     */
    static const bool is_complex = true;

    /**
     * For this data type, typedef the corresponding real type. Since this
     * specialization of the template is selected for number types
     * std::complex<T>, the real type is equal to the type used to store the
     * two components of the complex number.
     */
    typedef number real_type;

    /**
     * Return the complex-conjugate of the given number.
     */
    static
    std::complex<number> conjugate (const std::complex<number> &x);

    /**
     * Return the square of the absolute value of the given number. Since this
     * specialization of the general template is chosen for types equal to
     * std::complex, this function returns the product of a number and its
     * complex conjugate.
     */
    static
    real_type abs_square (const std::complex<number> &x);


    /**
     * Return the absolute value of a complex number.
     */
    static
    real_type abs (const std::complex<number> &x);
  };

  // --------------- inline and template functions ---------------- //

  inline bool is_nan (const double x)
  {
    return std::isnan(x);
  }

  inline bool is_finite (const double x)
  {
    return std::isfinite(x);
  }



  inline bool is_finite (const std::complex<double> &x)
  {
    // Check complex numbers for infinity
    // by testing real and imaginary part
    return ( is_finite (x.real())
             &&
             is_finite (x.imag()) );
  }



  inline bool is_finite (const std::complex<float> &x)
  {
    // Check complex numbers for infinity
    // by testing real and imaginary part
    return ( is_finite (x.real())
             &&
             is_finite (x.imag()) );
  }



  inline bool is_finite (const std::complex<long double> &x)
  {
    // Same for std::complex<long double>
    return ( is_finite (x.real())
             &&
             is_finite (x.imag()) );
  }


  template <typename number>
  DEAL_II_CUDA_HOST_DEV
  const number &
  NumberTraits<number>::conjugate (const number &x)
  {
    return x;
  }



  template <typename number>
  DEAL_II_CUDA_HOST_DEV
  typename NumberTraits<number>::real_type
  NumberTraits<number>::abs_square (const number &x)
  {
    return x * x;
  }



  template <typename number>
  typename NumberTraits<number>::real_type
  NumberTraits<number>::abs (const number &x)
  {
    return std::abs(x);
  }



  template <typename number>
  std::complex<number>
  NumberTraits<std::complex<number> >::conjugate (const std::complex<number> &x)
  {
    return std::conj(x);
  }



  template <typename number>
  typename NumberTraits<std::complex<number> >::real_type
  NumberTraits<std::complex<number> >::abs (const std::complex<number> &x)
  {
    return std::abs(x);
  }



  template <typename number>
  typename NumberTraits<std::complex<number> >::real_type
  NumberTraits<std::complex<number> >::abs_square (const std::complex<number> &x)
  {
    return std::norm (x);
  }

}

namespace internal
{
  /**
   * A test to see if it is possible to convert one number
   * type to the other.
   */
  template<typename From, typename To>
  struct is_explicitly_convertible
  {
    // Source: https://stackoverflow.com/a/16944130
  private:
    template<typename T>
    static void f(T);

    template<typename F, typename T>
    static constexpr auto test(int) ->
    decltype(f(static_cast<T>(std::declval<F>())),true)
    {
      return true;
    }

    template<typename F, typename T>
    static constexpr auto test(...) -> bool
    {
      return false;
    }

  public:

    static bool const value = test<From,To>(0);
  };

  /**
  * The structs below are needed since VectorizedArray<T1> is a POD-type without a constructor and
  * can be a template argument for SymmetricTensor<...,T2> where T2 would equal VectorizedArray<T1>.
  * Internally, in previous versions of deal.II, SymmetricTensor<...,T2> would make use of the constructor
  * of T2 leading to a compile-time error. However simply adding a constructor for VectorizedArray<T1>
  * breaks the POD-idioms needed elsewhere. Calls to constructors of T2 subsequently got replaced by a
  * call to internal::NumberType<T2> which then determines the right function to use by template deduction.
  * A detailled discussion can be found at https://github.com/dealii/dealii/pull/3967 . Also see
  * numbers.h for another specialization.
  */
  template <typename T>
  struct NumberType
  {
    static DEAL_II_CUDA_HOST_DEV const T &value (const T &t)
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
    template<typename F>
    static T
    value (const F &f,
           typename std::enable_if<
           !std::is_same<typename std::decay<T>::type,typename std::decay<F>::type>::value &&
           std::is_constructible<T,F>::value
           >::type * = nullptr)
    {
      return T(f);
    }

    // Type T is explicitly convertible (but not constructible) from F.
    template<typename F>
    static T
    value (const F &f,
           typename std::enable_if<
           !std::is_same<typename std::decay<T>::type,typename std::decay<F>::type>::value &&
           !std::is_constructible<T,F>::value &&
           is_explicitly_convertible<const F,T>::value
           >::type * = 0)
    {
      return static_cast<T>(f);
    }
  };

  template <typename T>
  struct NumberType<std::complex<T> >
  {
    static const std::complex<T> &value (const std::complex<T> &t)
    {
      return t;
    }

    static std::complex<T> value (const T &t)
    {
      return std::complex<T>(t);
    }

    // Facilitate cast from complex<double> to complex<float>
    template <typename U>
    static std::complex<T> value (const std::complex<U> &t)
    {
      return std::complex<T>(
               NumberType<T>::value(t.real()),
               NumberType<T>::value(t.imag()));
    }
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
