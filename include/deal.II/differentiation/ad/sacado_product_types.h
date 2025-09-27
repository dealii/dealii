// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_ad_sacado_product_types_h
#define dealii_differentiation_ad_sacado_product_types_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

#ifdef DEAL_II_TRILINOS_WITH_SACADO

#  include <Sacado.hpp>
// It appears that some versions of Trilinos do not directly or indirectly
// include all the headers for all forward and reverse Sacado AD types.
// So we directly include these both here as a precaution.
// Standard forward AD classes (templated)
#  include <Sacado_Fad_DFad.hpp>
// Reverse AD classes (templated)
#  include <Sacado_trad.hpp>

DEAL_II_NAMESPACE_OPEN


/* -------------- Sacado::Fad::DFad
 * (Differentiation::AD::NumberTypes::[sacado_fad/sacado_fad_fad])
 * -------------- */


namespace internal
{
  template <typename T>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, float>
  {
    using type = Sacado::Fad::DFad<T>;
  };

  template <typename T>
  struct ProductTypeImpl<float, Sacado::Fad::DFad<T>>
  {
    using type = Sacado::Fad::DFad<T>;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, double>
  {
    using type = Sacado::Fad::DFad<T>;
  };

  template <typename T>
  struct ProductTypeImpl<double, Sacado::Fad::DFad<T>>
  {
    using type = Sacado::Fad::DFad<T>;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, int>
  {
    using type = Sacado::Fad::DFad<T>;
  };

  template <typename T>
  struct ProductTypeImpl<int, Sacado::Fad::DFad<T>>
  {
    using type = Sacado::Fad::DFad<T>;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, Sacado::Fad::DFad<U>>
  {
    using type = Sacado::Fad::DFad<typename ProductType<T, U>::type>;
  };


  // Sacado::Fad::Dfad expression templates
  // We demote the result of the expression template operations
  // to a Sacado::Fad::Dfad itself. This is the only way to retain
  // consistency between the number type going into a complex chain of
  // (potentially branching) operations and that coming out of them.

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Fad::Expr<T>, U>
  {
    using type =
      typename ProductType<typename Sacado::Fad::Expr<T>::value_type, U>::type;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<T, Sacado::Fad::Expr<U>>
  {
    using type =
      typename ProductType<T, typename Sacado::Fad::Expr<U>::value_type>::type;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Fad::Expr<T>, Sacado::Fad::Expr<U>>
  {
    using type =
      typename ProductType<typename Sacado::Fad::Expr<T>::value_type,
                           typename Sacado::Fad::Expr<U>::value_type>::type;
  };

} // namespace internal


template <typename T>
struct EnableIfScalar<Sacado::Fad::DFad<T>>
{
  using type = Sacado::Fad::DFad<T>;
};

template <typename T>
struct EnableIfScalar<Sacado::Fad::Expr<T>>
{
  using type = typename Sacado::Fad::Expr<T>::value_type;
};


/* -------------- Sacado::Rad::ADvar
 * (Differentiation::AD::NumberTypes::[sacado_rad/sacado_rad_fad])
 * -------------- */


namespace internal
{
  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, float>
  {
    using type = Sacado::Rad::ADvar<T>;
  };

  template <typename T>
  struct ProductTypeImpl<float, Sacado::Rad::ADvar<T>>
  {
    using type = Sacado::Rad::ADvar<T>;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, double>
  {
    using type = Sacado::Rad::ADvar<T>;
  };

  template <typename T>
  struct ProductTypeImpl<double, Sacado::Rad::ADvar<T>>
  {
    using type = Sacado::Rad::ADvar<T>;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, int>
  {
    using type = Sacado::Rad::ADvar<T>;
  };

  template <typename T>
  struct ProductTypeImpl<int, Sacado::Rad::ADvar<T>>
  {
    using type = Sacado::Rad::ADvar<T>;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, Sacado::Rad::ADvar<U>>
  {
    using type = Sacado::Rad::ADvar<typename ProductType<T, U>::type>;
  };

  /* --- Sacado::Rad::ADvar: Temporary type --- */

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, float>
  {
    using type = Sacado::Rad::ADvari<T>;
  };

  template <typename T>
  struct ProductTypeImpl<float, Sacado::Rad::ADvari<T>>
  {
    using type = Sacado::Rad::ADvari<T>;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, double>
  {
    using type = Sacado::Rad::ADvari<T>;
  };

  template <typename T>
  struct ProductTypeImpl<double, Sacado::Rad::ADvari<T>>
  {
    using type = Sacado::Rad::ADvari<T>;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, int>
  {
    using type = Sacado::Rad::ADvari<T>;
  };

  template <typename T>
  struct ProductTypeImpl<int, Sacado::Rad::ADvari<T>>
  {
    using type = Sacado::Rad::ADvari<T>;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, Sacado::Rad::ADvari<U>>
  {
    using type = Sacado::Rad::ADvari<typename ProductType<T, U>::type>;
  };

  /* --- Sacado::Rad::ADvar: Main and temporary type --- */

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, Sacado::Rad::ADvari<U>>
  {
    using type = Sacado::Rad::ADvar<typename ProductType<T, U>::type>;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, Sacado::Rad::ADvar<U>>
  {
    using type = Sacado::Rad::ADvar<typename ProductType<T, U>::type>;
  };

} // namespace internal


template <typename T>
struct EnableIfScalar<Sacado::Rad::ADvar<T>>
{
  using type = Sacado::Rad::ADvar<T>;
};


template <typename T>
struct EnableIfScalar<Sacado::Rad::ADvari<T>>
{
  using type = Sacado::Rad::ADvari<T>;
};


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_SACADO

#endif
