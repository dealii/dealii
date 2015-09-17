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


/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of a
 * sacado type with a different non-sacado type.
 *
 * @relates EnableIfScalar
 * @relates ProductType
 */

template <typename T, typename U>
typename ProductType<Sacado::Fad::DFad<T>, typename EnableIfScalar<U>::type>::type
inline
operator*(const Sacado::Fad::DFad<T> &left, const U &right)
{
  typedef typename ProductType<Sacado::Fad::DFad<T>, U>::type result_type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}


/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of non-sacado type with a sacado type.
 *
 * @relates EnableIfScalar
 * @relates ProductType
 */

template <typename T, typename U>
typename ProductType<typename EnableIfScalar<T>::type, Sacado::Fad::DFad<U> >::type
inline
operator*(const T &left, const Sacado::Fad::DFad<U> &right)
{
  typedef typename ProductType<T, Sacado::Fad::DFad<U> >::type result_type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * Provide an <tt>operator*</tt> for a scalar multiplication of mixed
 * sacado types.
 *
 * @relates EnableIfScalar
 * @relates ProductType
 */

template <typename T, typename U>
typename ProductType<Sacado::Fad::DFad<T>, Sacado::Fad::DFad<U> >::type
inline
operator*(const Sacado::Fad::DFad<T> &left, const Sacado::Fad::DFad<U> &right)
{
  typedef typename ProductType<Sacado::Fad::DFad<T>, Sacado::Fad::DFad<U> >::type result_type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

#endif // DEAL_II_WITH_TRILINOS

DEAL_II_NAMESPACE_CLOSE

#endif
