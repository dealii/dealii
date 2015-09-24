// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

#ifndef dealii__sacado_product_type_h
#define dealii__sacado_product_type_h

#include <deal.II/base/config.h>
#include <deal.II/base/template_constraints.h>

#ifdef DEAL_II_WITH_TRILINOS
#include "Sacado.hpp"

DEAL_II_NAMESPACE_OPEN

template <typename T>
struct ProductType<Sacado::Fad::DFad<T>, float>
{
  typedef Sacado::Fad::DFad<T> type;
};

template <typename T>
struct ProductType<float, Sacado::Fad::DFad<T> >
{
  typedef Sacado::Fad::DFad<T> type;
};

template <typename T>
struct ProductType<Sacado::Fad::DFad<T>, double>
{
  typedef Sacado::Fad::DFad<T> type;
};

template <typename T>
struct ProductType<double, Sacado::Fad::DFad<T> >
{
  typedef Sacado::Fad::DFad<T> type;
};

template <typename T>
struct ProductType<Sacado::Fad::DFad<T>, int>
{
  typedef Sacado::Fad::DFad<T> type;
};

template <typename T>
struct ProductType<int, Sacado::Fad::DFad<T> >
{
  typedef Sacado::Fad::DFad<T> type;
};

template <typename T, typename U>
struct ProductType<Sacado::Fad::DFad<T>, Sacado::Fad::DFad<U> >
{
  typedef Sacado::Fad::DFad<typename ProductType<T,U>::type > type;
};

template <typename T>
struct EnableIfScalar<Sacado::Fad::DFad<T> >
{
  typedef Sacado::Fad::DFad<T> type;
};

#endif // DEAL_II_WITH_TRILINOS

DEAL_II_NAMESPACE_CLOSE

#endif
