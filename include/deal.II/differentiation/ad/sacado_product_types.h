// ---------------------------------------------------------------------
//
// Copyright (C) 2015, 2017 by the deal.II authors
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

#ifndef dealii_differentiation_ad_sacado_product_types_h
#define dealii_differentiation_ad_sacado_product_types_h

#include <deal.II/base/config.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>

#ifdef DEAL_II_WITH_TRILINOS

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <Sacado.hpp>
// It appears that some versions of Trilinos do not directly or indirectly
// include all the headers for all forward and reverse Sacado AD types.
// So we directly include these both here as a precaution.
// Standard forward AD classes (templated)
#include <Sacado_Fad_DFad.hpp>
// Reverse AD classes (templated)
#include <Sacado_trad.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN


/* -------------- Sacado::Fad::DFad (Differentiation::AD::NumberTypes::[sacado_fad/sacado_fad_fad]) -------------- */


namespace internal
{

  template <typename T>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, float>
  {
    typedef Sacado::Fad::DFad<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<float, Sacado::Fad::DFad<T> >
  {
    typedef Sacado::Fad::DFad<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, double>
  {
    typedef Sacado::Fad::DFad<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<double, Sacado::Fad::DFad<T> >
  {
    typedef Sacado::Fad::DFad<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, int>
  {
    typedef Sacado::Fad::DFad<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<int, Sacado::Fad::DFad<T> >
  {
    typedef Sacado::Fad::DFad<T> type;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Fad::DFad<T>, Sacado::Fad::DFad<U> >
  {
    typedef Sacado::Fad::DFad<typename ProductType<T,U>::type> type;
  };

}


template <typename T>
struct EnableIfScalar<Sacado::Fad::DFad<T> >
{
  typedef Sacado::Fad::DFad<T> type;
};


/* -------------- Sacado::Rad::ADvar (Differentiation::AD::NumberTypes::[sacado_rad/sacado_rad_fad]) -------------- */


namespace internal
{

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, float>
  {
    typedef Sacado::Rad::ADvar<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<float, Sacado::Rad::ADvar<T> >
  {
    typedef Sacado::Rad::ADvar<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, double>
  {
    typedef Sacado::Rad::ADvar<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<double, Sacado::Rad::ADvar<T> >
  {
    typedef Sacado::Rad::ADvar<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, int>
  {
    typedef Sacado::Rad::ADvar<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<int, Sacado::Rad::ADvar<T> >
  {
    typedef Sacado::Rad::ADvar<T> type;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, Sacado::Rad::ADvar<U> >
  {
    typedef Sacado::Rad::ADvar<typename ProductType<T,U>::type> type;
  };

  /* --- Sacado::Rad::ADvar: Temporary type --- */

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, float>
  {
    typedef Sacado::Rad::ADvari<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<float, Sacado::Rad::ADvari<T> >
  {
    typedef Sacado::Rad::ADvari<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, double>
  {
    typedef Sacado::Rad::ADvari<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<double, Sacado::Rad::ADvari<T> >
  {
    typedef Sacado::Rad::ADvari<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, int>
  {
    typedef Sacado::Rad::ADvari<T> type;
  };

  template <typename T>
  struct ProductTypeImpl<int, Sacado::Rad::ADvari<T> >
  {
    typedef Sacado::Rad::ADvari<T> type;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, Sacado::Rad::ADvari<U> >
  {
    typedef Sacado::Rad::ADvari<typename ProductType<T,U>::type> type;
  };

  /* --- Sacado::Rad::ADvar: Main and temporary type --- */

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvar<T>, Sacado::Rad::ADvari<U> >
  {
    typedef Sacado::Rad::ADvar<typename ProductType<T,U>::type> type;
  };

  template <typename T, typename U>
  struct ProductTypeImpl<Sacado::Rad::ADvari<T>, Sacado::Rad::ADvar<U> >
  {
    typedef Sacado::Rad::ADvar<typename ProductType<T,U>::type> type;
  };

}


template <typename T>
struct EnableIfScalar<Sacado::Rad::ADvar<T> >
{
  typedef Sacado::Rad::ADvar<T> type;
};


template <typename T>
struct EnableIfScalar<Sacado::Rad::ADvari<T> >
{
  typedef Sacado::Rad::ADvari<T> type;
};


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

#endif
