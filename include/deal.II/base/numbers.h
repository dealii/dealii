// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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

#ifndef __deal2__numbers_h
#define __deal2__numbers_h


#include <deal.II/base/config.h>
#include <deal.II/base/types.h>
#include <complex>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for the declaration of universal constants. Since the
 * availability in <tt>math.h</tt> is not always guaranteed, we put
 * them here. Since this file is included by <tt>base/config.h</tt>,
 * they are available to the whole library.
 *
 * The constants defined here are a subset of the <tt>M_XXX</tt> constants
 * sometimes declared in the system include file <tt>math.h</tt>, but without
 * the prefix <tt>M_</tt>.
 *
 * In addition to that, we declare  <tt>invalid_unsigned_int</tt> to be the
 * largest unsigned integer representable; this value is widely used in
 * the library as a marker for an invalid index, an invalid size of an
 * array, and similar purposes.
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
   * Return @p true if the given
   * value is a finite floating
   * point number, i.e. is neither
   * plus or minus infinity nor NaN
   * (not a number).
   *
   * Note that the argument type of
   * this function is
   * <code>double</code>. In other
   * words, if you give a very
   * large number of type
   * <code>long double</code>, this
   * function may return
   * <code>false</code> even if the
   * number is finite with respect
   * to type <code>long
   * double</code>.
   */
  bool is_finite (const double x);

  /**
   * Return @p true if real and
   * imaginary parts of the given
   * complex number are finite.
   */
  bool is_finite (const std::complex<double> &x);

  /**
   * Return @p true if real and
   * imaginary parts of the given
   * complex number are finite.
   */
  bool is_finite (const std::complex<float> &x);

  /**
   * Return @p true if real and
   * imaginary parts of the given
   * complex number are finite.
   *
   * Again may not work correctly if
   * real or imaginary parts are very
   * large numbers that are infinite in
   * terms of <code>double</code>, but
   * finite with respect to
   * <code>long double</code>.
   */
  bool is_finite (const std::complex<long double> &x);

  /**
   * A structure that, together with
   * its partial specializations
   * NumberTraits<std::complex<number> >,
   * provides traits and member
   * functions that make it possible
   * to write templates that work on
   * both real number types and
   * complex number types. This
   * template is mostly used to
   * implement linear algebra classes
   * such as vectors and matrices
   * that work for both real and
   * complex numbers.
   *
   * @author Wolfgang Bangerth, 2007
   */
  template <typename number>
  struct NumberTraits
  {
    /**
     * A flag that specifies whether the
     * template type given to this class is
     * complex or real. Since the general
     * template is selected for non-complex
     * types, the answer is
     * <code>false</code>.
     */
    static const bool is_complex = false;

    /**
     * For this data type, typedef the
     * corresponding real type. Since the
     * general template is selected for all
     * data types that are not
     * specializations of std::complex<T>,
     * the underlying type must be
     * real-values, so the real_type is
     * equal to the underlying type.
     */
    typedef number real_type;

    /**
     * Return the complex-conjugate of the
     * given number. Since the general
     * template is selected if number is
     * not a complex data type, this
     * function simply returns the given
     * number.
     */
    static
    const number &conjugate (const number &x);

    /**
     * Return the square of the absolute
     * value of the given number. Since the
     * general template is chosen for types
     * not equal to std::complex, this
     * function simply returns the square
     * of the given number.
     */
    static
    real_type abs_square (const number &x);

    /**
     * Return the absolute value of a
     * number.
     */
    static
    real_type abs (const number &x);
  };


  /**
   * Specialization of the general
   * NumberTraits class that provides the
   * relevant information if the underlying
   * data type is std::complex<T>.
   *
   * @author Wolfgang Bangerth, 2007
   */
  template <typename number>
  struct NumberTraits<std::complex<number> >
  {
    /**
     * A flag that specifies whether the
     * template type given to this class is
     * complex or real. Since this
     * specialization of the general
     * template is selected for complex
     * types, the answer is
     * <code>true</code>.
     */
    static const bool is_complex = true;

    /**
     * For this data type, typedef the
     * corresponding real type. Since this
     * specialization of the template is
     * selected for number types
     * std::complex<T>, the real type is
     * equal to the type used to store the
     * two components of the complex
     * number.
     */
    typedef number real_type;

    /**
     * Return the complex-conjugate of the
     * given number.
     */
    static
    std::complex<number> conjugate (const std::complex<number> &x);

    /**
     * Return the square of the absolute
     * value of the given number. Since
     * this specialization of the general
     * template is chosen for types equal
     * to std::complex, this function
     * returns the product of a number and
     * its complex conjugate.
     */
    static
    real_type abs_square (const std::complex<number> &x);


    /**
     * Return the absolute value of a
     * complex number.
     */
    static
    real_type abs (const std::complex<number> &x);
  };

}


//TODO[WB]: eventually remove this namespace alias
/*
 * Namespace alias with the old name for the numbers namespace. The namespace
 * was originally called numbers, but has since been renamed to
 * dealii::numbers when everything was moved into namespace dealii.
 *
 * @deprecated
 */
namespace deal_II_numbers = numbers;


// --------------- inline and template functions ---------------- //

namespace numbers
{
  template <typename number>
  const number &
  NumberTraits<number>::conjugate (const number &x)
  {
    return x;
  }



  template <typename number>
  typename NumberTraits<number>::real_type
  NumberTraits<number>::abs_square (const number &x)
  {
    return x * x;
  }



  template <typename number>
  typename NumberTraits<number>::real_type
  NumberTraits<number>::abs (const number &x)
  {
    return std::fabs(x);
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



DEAL_II_NAMESPACE_CLOSE

#endif
