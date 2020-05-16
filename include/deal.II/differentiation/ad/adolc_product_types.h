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

#ifndef dealii_differentiation_ad_adolc_product_types_h
#define dealii_differentiation_ad_adolc_product_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ADOLC

#  include <deal.II/base/template_constraints.h>

#  include <adolc/adouble.h> // Taped double
#  include <adolc/adtl.h>    // Tapeless double

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

#endif // DEAL_II_WITH_ADOLC

#endif
