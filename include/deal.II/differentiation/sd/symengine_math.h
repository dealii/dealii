// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_sd_symengine_math_h
#define dealii_differentiation_sd_symengine_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_number_types.h>

#  include <type_traits>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    // Mathematical functions:
    //
    // It is necessary that all computed inputs to SymEngine functions are
    // of type SymEngine::RCP<SymEngine::Basic>. So we instead simply offer a
    // unified interface to Expression types and with other permutations of
    // numbers we convert between them.
    //
    // For a full list of functions that we (ultimately) expect to support, see
    // http://www.cplusplus.com/reference/cmath/. Those that are currently
    // supported are extracted from symengine/{functions,pow}.h;
    // symengine/type_codes.inc.


    /**
     * @name Power functions
     */
    /** @{ */

    /**
     * Return a symbolic number that represents a @p base value raised to the power of
     * an @p exponent.
     *
     * Mimics the function <code> std::%pow(base,exponent) </code> using the
     * standard math library.
     */
    Expression
    pow(const Expression &base, const Expression &exponent);

    /**
     * Return a symbolic number that represents a @p base value raised to the power of
     * an @p exponent.
     *
     * Mimics the function <code> std::%pow(base,exponent) </code> using the
     * standard math library.
     *
     * This variant is used when the @p exponent is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    pow(const Expression &base, const NumberType &exponent)
    {
      // Call other implementation
      return pow(base, Expression(exponent));
    }

    /**
     * Return a symbolic number that represents a @p base value raised to the power of
     * an @p exponent.
     *
     * Mimics the function <code> std::%pow(base,exponent) </code> using the
     * standard math library.
     *
     * This variant is used when the @p base is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    pow(const NumberType &base, const Expression &exponent)
    {
      // Call other implementation
      return pow(Expression(base), exponent);
    }

    /**
     * Return a symbolic number that represents the square root of some value @p x.
     *
     * Mimics the function <code> std::%sqrt(x) </code> using the standard math
     * library.
     */
    Expression
    sqrt(const Expression &x);

    /**
     * Return a symbolic number that represents the cubic root of some value @p x.
     *
     * Mimics the function <code> std::%cbrt(x) </code> using the standard math
     * library.
     */
    Expression
    cbrt(const Expression &x);

    /**
     * Return a symbolic number that represents the Euler constant
     * $e \approx 2.71828$ raised to the given @p exponent.
     *
     * Mimics the function <code> std::%exp(exponent) </code> using the standard
     * math library.
     */
    Expression
    exp(const Expression &exponent);

    /**
     * Return a symbolic number that represents the natural logarithm of a value @p x.
     *
     * Mimics the function <code> std::%log(x) </code> using the standard math
     * library.
     */
    Expression
    log(const Expression &x);

    /**
     * Return a symbolic number that represents the logarithm of a value @p x taken with
     * respect to a @p base number.
     *
     * Mimics the function <code> std::%log(x,base) </code> using the standard
     * math library.
     */
    Expression
    log(const Expression &x, const Expression &base);

    /**
     * Return a symbolic number that represents the logarithm of a value @p x taken with
     * respect to a @p base number.
     *
     * Mimics the function <code> std::%log(x,base) </code> using the standard
     * math library.
     *
     * This variant is used when the @p base is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    log(const Expression &x, const NumberType &base)
    {
      // Call other implementation
      return log(x, Expression(base));
    }

    /**
     * Return a symbolic number that represents the logarithm of a value @p x taken with
     * respect to a @p base number.
     *
     * Mimics the function <code> std::%log(x,base) </code> using the standard
     * math library.
     *
     * This variant is used when the @p value is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    log(const NumberType &x, const Expression &base)
    {
      // Call other implementation
      return log(Expression(x), base);
    }

    /**
     * Return a symbolic number that represents the base 10 logarithm of a value @p x.
     *
     * Mimics the function <code> std::%log10(x) </code> using the standard math
     * library.
     */
    Expression
    log10(const Expression &x);

    /** @} */

    /**
     * @name Trigonometric functions
     */
    /** @{ */

    /**
     * Return a symbolic number that represents the sine function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%sin(x) </code> using the standard math
     * library.
     */
    Expression
    sin(const Expression &x);

    /**
     * Return a symbolic number that represents the cosine function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%cos(x) </code> using the standard math
     * library.
     */
    Expression
    cos(const Expression &x);

    /**
     * Return a symbolic number that represents the tangent function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%tan(x) </code> using the standard math
     * library.
     */
    Expression
    tan(const Expression &x);

    /**
     * Return a symbolic number that represents the cosecant function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%sin(x) </code> using the standard
     * math library.
     */
    Expression
    csc(const Expression &x);

    /**
     * Return a symbolic number that represents the secant function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%cos(x) </code> using the standard
     * math library.
     */
    Expression
    sec(const Expression &x);

    /**
     * Return a symbolic number that represents the cotangent function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%tan(x) </code> using the standard
     * math library.
     */
    Expression
    cot(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse sine function with
     * the
     * given argument @p x.
     *
     * Mimics the function <code> std::%asin(x) </code> using the standard math
     * library.
     */
    Expression
    asin(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse cosine function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%acos(x) </code> using the standard math
     * library.
     */
    Expression
    acos(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse tangent function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%atan(x) </code> using the standard math
     * library.
     */
    Expression
    atan(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse tangent function
     * with the
     * given arguments @p x and @p y.
     *
     * Mimics the function <code> std::%atan2(y,x) </code> using the standard
     * math library.
     */
    Expression
    atan2(const Expression &y, const Expression &x);

    /**
     * Return a symbolic number that represents the inverse tangent function
     * with the
     * given arguments @p x and @p y.
     *
     * Mimics the function <code> std::%atan2(y,x) </code> using the standard
     * math library.
     *
     * This variant is used when the numerator @p y is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    atan2(const NumberType &y, const Expression &x)
    {
      // Call other implementation
      return atan2(Expression(y), x);
    }

    /**
     * Return a symbolic number that represents the inverse tangent function
     * with the
     * given arguments @p x and @p y.
     *
     * Mimics the function <code> std::%atan2(y,x) </code> using the standard
     * math library.
     *
     * This variant is used when the denominator @p x is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    atan2(const Expression &y, const NumberType &x)
    {
      // Call other implementation
      return atan2(y, Expression(x));
    }

    /**
     * Return a symbolic number that represents the inverse cosecant function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%asin(x) </code> using the standard
     * math library.
     */
    Expression
    acsc(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse secant function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%acos(x) </code> using the standard
     * math library.
     */
    Expression
    asec(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse cotangent function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%atan(x) </code> using the standard
     * math library.
     */
    Expression
    acot(const Expression &x);

    /** @} */

    /**
     * @name Hyperbolic trigonometric functions
     */
    /** @{ */

    /**
     * Return a symbolic number that represents the hyperbolic sine function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%sinh(x) </code> using the standard math
     * library.
     */
    Expression
    sinh(const Expression &x);

    /**
     * Return a symbolic number that represents the hyperbolic cosine function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%cosh(x) </code> using the standard math
     * library.
     */
    Expression
    cosh(const Expression &x);

    /**
     * Return a symbolic number that represents the hyperbolic tangent function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%tanh(x) </code> using the standard math
     * library.
     */
    Expression
    tanh(const Expression &x);

    /**
     * Return a symbolic number that represents the hyperbolic cosecant
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%sinh(x) </code> using the standard
     * math library.
     */
    Expression
    csch(const Expression &x);

    /**
     * Return a symbolic number that represents the hyperbolic secant function
     * with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%cosh(x) </code> using the standard
     * math library.
     */
    Expression
    sech(const Expression &x);

    /**
     * Return a symbolic number that represents the hyperbolic cotangent
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%tanh(x) </code> using the standard
     * math library.
     */
    Expression
    coth(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse hyperbolic sine
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%asinh(x) </code> using the standard math
     * library.
     */
    Expression
    asinh(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse hyperbolic cosine
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%acosh(x) </code> using the standard math
     * library.
     */
    Expression
    acosh(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse hyperbolic tangent
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%atanh(x) </code> using the standard math
     * library.
     */
    Expression
    atanh(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse hyperbolic cosecant
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%asin(x) </code> using the standard
     * math library.
     */
    Expression
    acsch(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse hyperbolic secant
     * function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%acos(x) </code> using the standard
     * math library.
     */
    Expression
    asech(const Expression &x);

    /**
     * Return a symbolic number that represents the inverse hyperbolic
     * cotangent function with the
     * given argument @p x.
     *
     * Mimics the function <code> 1.0/std::%atan(x) </code> using the standard
     * math library.
     */
    Expression
    acoth(const Expression &x);

    /** @} */

    /**
     * @name Other functions
     */
    /** @{ */

    /**
     * Return a symbolic number that represents the absolute value of value @p x.
     *
     * Mimics the function <code> std::%abs(x) </code> using the standard math
     * library.
     */
    Expression
    abs(const Expression &x);


    /**
     * Return a symbolic number that represents the absolute value of value @p x.
     *
     * Mimics the function <code> std::%fabs(x) </code> using the standard math
     * library.
     */
    Expression
    fabs(const Expression &x);


    /**
     * Return a symbolic number that represents the sign of value @p x.
     *
     * Although there is no such function in the standard library, it mimics
     * the function <code> boost::sign(x) </code> using the
     * boost math library.
     */
    Expression
    sign(const Expression &x);


    /**
     * Return a symbolic number that represents the @p value of the first
     * argument that takes the @p sign of the second argument.
     *
     * Mimics the function <code> std::%copysign(value, sign) </code> using
     * the standard math library.
     */
    Expression
    copysign(const Expression &value, const Expression &sign);


    /**
     * Return a symbolic number that represents the floor of value @p x.
     *
     * Mimics the function <code> std::%floor(x) </code> using the standard math
     * library.
     */
    Expression
    floor(const Expression &x);


    /**
     * Return a symbolic number that represents the ceiling of value @p x.
     *
     * Mimics the function <code> std::%ceil(x) </code> using the standard math
     * library.
     */
    Expression
    ceil(const Expression &x);

    /**
     * Return a symbolic number that represents the maximum of two
     * values @p a and @p b.
     *
     * Mimics the function <code> std::%max(a,b) </code> using the standard math
     * library.
     */
    Expression
    max(const Expression &a, const Expression &b);

    /**
     * Return a symbolic number that represents the maximum of two
     * values @p a and @p b.
     *
     * Mimics the function <code> std::%max(a,b) </code> using the standard math
     * library.
     *
     * This variant is used when @p b is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    max(const Expression &a, const NumberType &b)
    {
      // Call other implementation
      return max(a, Expression(b));
    }

    /**
     * Return a symbolic number that represents the maximum of two
     * values @p a and @p b.
     *
     * Mimics the function <code> std::%max(a,b) </code> using the standard math
     * library.
     *
     * This variant is used when @p a is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    max(const NumberType &a, const Expression &b)
    {
      // Call other implementation
      return max(Expression(a), b);
    }

    /**
     * Return a symbolic number that represents the minimum of two
     * values @p a and @p b.
     *
     * Mimics the function <code> std::%min(a,b) </code> using the standard math
     * library.
     */
    Expression
    min(const Expression &a, const Expression &b);

    /**
     * Return a symbolic number that represents the minimum of two
     * values @p a and @p b.
     *
     * Mimics the function <code> std::%min(a,b) </code> using the standard math
     * library.
     *
     * This variant is used when @p b is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    min(const Expression &a, const NumberType &b)
    {
      // Call other implementation
      return min(a, Expression(b));
    }

    /**
     * Return a symbolic number that represents the minimum of two
     * values @p a and @p b.
     *
     * Mimics the function <code> std::%min(a,b) </code> using the standard math
     * library.
     *
     * This variant is used when @p a is not a Expression.
     */
    template <
      typename NumberType,
      typename = std::enable_if_t<!std::is_same_v<NumberType, Expression>>>
    Expression
    min(const NumberType &a, const Expression &b)
    {
      // Call other implementation
      return min(Expression(a), b);
    }

    /**
     * Return a symbolic number that represents error function with the
     * given argument @p x.
     *
     * Mimics the function <code> std::%erf(x) </code> using the standard math
     * library.
     */
    Expression
    erf(const Expression &x);

    /**
     * Return a symbolic number that represents complimentary error function
     * with the given argument @p x.
     *
     * Mimics the function <code> std::%erfc(x) </code> using the standard math
     * library.
     */
    Expression
    erfc(const Expression &x);

    /** @} */

  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE


#  ifndef DOXYGEN

// Import math operations into standard namespace.
// This gives us the ability to use them within the Tensor class.
namespace std
{
  /**
   * Expose SymEngine wrapper math functions
   */

#    define DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(func) \
      using dealii::Differentiation::SD::func;

#    define DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION(func) \
      using dealii::Differentiation::SD::func;

  DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION(pow)
  DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION(max)
  DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION(min)
  DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION(copysign)

  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(exp)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(log10)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(sqrt)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(cbrt)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(erf)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(erfc)
  // Note: Both unary and binary versions
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(log)

  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(abs)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(fabs)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(sign)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(floor)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(ceil)

  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(sin)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(cos)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(tan)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(csc)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(sec)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(cot)

  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(asin)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(acos)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(atan)
  DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION(atan2)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(acsc)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(asec)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(acot)

  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(sinh)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(cosh)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(tanh)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(csch)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(sech)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(coth)

  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(asinh)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(acosh)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(atanh)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(acsch)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(asech)
  DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION(acoth)

#    undef DEAL_II_EXPOSE_SYMENGINE_BINARY_MATH_FUNCTION
#    undef DEAL_II_EXPOSE_SYMENGINE_UNARY_MATH_FUNCTION

} // namespace std

#  endif // DOXYGEN

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
