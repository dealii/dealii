// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_adolc_product_types_h
#define dealii_differentiation_ad_adolc_product_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <deal.II/base/template_constraints.h>

#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double

#  include <complex>


DEAL_II_NAMESPACE_OPEN


/* ------ ADOL-C taped (Differentiation::AD::NumberTypes::adolc_taped) ----- */


namespace internal
{
  template <>
  struct ProductTypeImpl<adouble, adouble>
  {
    using type = adouble;
  };

  // Typedefs for "adub"s are all that's necessary to ensure that no temporary
  // ADOL-C types "adub" are created when a scalar product is performed. If this
  // is not done, then intermediate tensors are filled with unconstructable
  // types.
  template <>
  struct ProductTypeImpl<adub, adouble>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<adouble, adub>
  {
    using type = adouble;
  };

  /* --- Double --- */

  template <>
  struct ProductTypeImpl<double, adouble>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<adouble, double>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<double, adub>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<adub, double>
  {
    using type = adouble;
  };

  /* --- Float --- */

  template <>
  struct ProductTypeImpl<float, adouble>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<adouble, float>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<float, adub>
  {
    using type = adouble;
  };

  template <>
  struct ProductTypeImpl<adub, float>
  {
    using type = adouble;
  };

  /* --- Complex double --- */

  template <>
  struct ProductTypeImpl<std::complex<double>, std::complex<adouble>>
  {
    using type = std::complex<adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adouble>, std::complex<double>>
  {
    using type = std::complex<adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adouble>, std::complex<adouble>>
  {
    using type = std::complex<adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adub>, std::complex<adouble>>
  {
    using type = std::complex<adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adouble>, std::complex<adub>>
  {
    using type = std::complex<adouble>;
  };

  /* --- Complex float --- */

  template <>
  struct ProductTypeImpl<std::complex<float>, std::complex<adouble>>
  {
    using type = std::complex<adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adouble>, std::complex<float>>
  {
    using type = std::complex<adouble>;
  };

} // namespace internal

template <>
struct EnableIfScalar<adouble>
{
  using type = adouble;
};

template <>
struct EnableIfScalar<std::complex<adouble>>
{
  using type = std::complex<adouble>;
};


template <>
struct EnableIfScalar<adub>
{
  using type = adouble;
};


template <>
struct EnableIfScalar<std::complex<adub>>
{
  using type = std::complex<adouble>;
};


/* -- ADOL-C tapeless (Differentiation::AD::NumberTypes::adolc_tapeless) -- */


namespace internal
{
  /* --- Double --- */

  template <>
  struct ProductTypeImpl<double, adtl::adouble>
  {
    using type = adtl::adouble;
  };

  template <>
  struct ProductTypeImpl<adtl::adouble, double>
  {
    using type = adtl::adouble;
  };

  template <>
  struct ProductTypeImpl<adtl::adouble, adtl::adouble>
  {
    using type = adtl::adouble;
  };

  /* --- Float --- */

  template <>
  struct ProductTypeImpl<float, adtl::adouble>
  {
    using type = adtl::adouble;
  };

  template <>
  struct ProductTypeImpl<adtl::adouble, float>
  {
    using type = adtl::adouble;
  };

  /* --- Complex double --- */

  template <>
  struct ProductTypeImpl<std::complex<double>, std::complex<adtl::adouble>>
  {
    using type = std::complex<adtl::adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adtl::adouble>, std::complex<double>>
  {
    using type = std::complex<adtl::adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adtl::adouble>,
                         std::complex<adtl::adouble>>
  {
    using type = std::complex<adtl::adouble>;
  };

  /* --- Complex float --- */

  template <>
  struct ProductTypeImpl<std::complex<float>, std::complex<adtl::adouble>>
  {
    using type = std::complex<adtl::adouble>;
  };

  template <>
  struct ProductTypeImpl<std::complex<adtl::adouble>, std::complex<float>>
  {
    using type = std::complex<adtl::adouble>;
  };

} // namespace internal


template <>
struct EnableIfScalar<adtl::adouble>
{
  using type = adtl::adouble;
};


template <>
struct EnableIfScalar<std::complex<adtl::adouble>>
{
  using type = std::complex<adtl::adouble>;
};


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_ADOLC

#endif
