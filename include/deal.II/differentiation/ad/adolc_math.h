// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#ifndef dealii_differentiation_ad_adolc_math_h
#define dealii_differentiation_ad_adolc_math_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double
#  include <adolc/internal/adolc_settings.h>
#  include <adolc/internal/adubfunc.h> // Taped double math functions


#  ifndef DOXYGEN

/**
 * Import ADOL-C math operations into standard namespace. This gives us the
 * ability to use them within the Tensor class, and it also allows the user
 * to write generic code and switch between AD number types.
 *
 * The math functions to be exposed come from the ADOL-C's taped
 * (anonymous) and tapeless namespaces.
 */
namespace std
{
  /**
   * Make a unary function with one name available under a different name
   */
#    define DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION_COPY(func_to, func_from) \
      inline adouble func_to(const adouble &x)                                \
      {                                                                       \
        return func_from(static_cast<const badouble &>(x));                   \
      }                                                                       \
      inline adtl::adouble func_to(const adtl::adouble &x)                    \
      {                                                                       \
        return adtl::func_from(x);                                            \
      }

  /**
   * Expose a unary function
   */
#    define DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(func) \
      DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION_COPY(func, func)

  /**
   * Make a binary function with one name available under a different name
   */
#    define DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_COPY(func_to, func_from) \
      inline adouble func_to(const adouble &x, const adouble &y)               \
      {                                                                        \
        return func_from(static_cast<const badouble &>(x),                     \
                         static_cast<const badouble &>(y));                    \
      }                                                                        \
      inline adouble func_to(const double x, const adouble &y)                 \
      {                                                                        \
        return func_from(x, static_cast<const badouble &>(y));                 \
      }                                                                        \
      inline adouble func_to(const adouble &x, const double y)                 \
      {                                                                        \
        return func_from(static_cast<const badouble &>(x), y);                 \
      }                                                                        \
      inline adtl::adouble func_to(const adtl::adouble &x,                     \
                                   const adtl::adouble &y)                     \
      {                                                                        \
        return adtl::func_from(x, y);                                          \
      }                                                                        \
      inline adtl::adouble func_to(const double x, const adtl::adouble &y)     \
      {                                                                        \
        return adtl::func_from(x, y);                                          \
      }                                                                        \
      inline adtl::adouble func_to(const adtl::adouble &x, const double y)     \
      {                                                                        \
        return adtl::func_from(x, y);                                          \
      }

  /**
   * Expose a binary function
   */
#    define DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION(func) \
      DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_COPY(func, func)

  /**
   * Expose a binary function
   */
#    define DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_2(func) \
      inline adouble func(const adouble &x, const adouble &y) \
      {                                                       \
        return func(static_cast<const badouble &>(x),         \
                    static_cast<const badouble &>(y));        \
      }                                                       \
      inline adtl::adouble func(const adtl::adouble &x,       \
                                const adtl::adouble &y)       \
      {                                                       \
        return adtl::func(x, y);                              \
      }

  // See
  // https://gitlab.com/adol-c/adol-c/blob/master/ADOL-C/include/adolc/adouble.h
  // https://gitlab.com/adol-c/adol-c/blob/master/ADOL-C/include/adolc/internal/paramfunc.h

  DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION(pow)
  DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION(fmax)
  DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_COPY(max, fmax)
  DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION(fmin)
  DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_COPY(min, fmin)

  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(exp)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(log)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(log10)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(sqrt)
#    if defined(DEAL_II_ADOLC_WITH_ATRIG_ERF)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(erf)
  inline adouble
  erfc(const adouble &x)
  {
    return 1.0 - std::erf(x);
  }
  inline adtl::adouble
  erfc(const adtl::adouble &x)
  {
    return 1.0 - std::erf(x);
  }
#    endif

  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(fabs)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION_COPY(abs, fabs)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(ceil)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(floor)

  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(sin)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(cos)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(tan)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(asin)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(acos)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(atan)
  DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_2(atan2)

  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(sinh)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(cosh)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(tanh)
#    if defined(DEAL_II_ADOLC_WITH_ATRIG_ERF)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(asinh)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(acosh)
  DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION(atanh)
#    endif

#    undef DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_2
#    undef DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION
#    undef DEAL_II_EXPOSE_ADOLC_BINARY_MATH_FUNCTION_COPY
#    undef DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION
#    undef DEAL_II_EXPOSE_ADOLC_UNARY_MATH_FUNCTION_COPY

} // namespace std

#  endif // DOXYGEN

#endif // DEAL_II_WITH_ADOLC

#endif
