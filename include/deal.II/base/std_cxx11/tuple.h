// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#ifndef __deal2__std_cxx11_tuple_h
#define __deal2__std_cxx11_tuple_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#  include <tuple>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using std::tuple;
  using std::get;
  using std::tuple_size;
  using std::tuple_element;
}
DEAL_II_NAMESPACE_CLOSE

#else

#include <boost/tuple/tuple.hpp>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using boost::tuple;
  using boost::get;

				   // boost::tuples::length has been renamed
				   // by the standard to std::tuple_size
  template <typename T>
  struct tuple_size 
  {
      static const std::size_t value = boost::tuples::length<T>::value;
  };

				   // similarly, boost::tuples::element has
				   // been renamed by the standard to
				   // std::tuple_element
  template <int N, typename T>
  struct tuple_element
  {
      typedef typename boost::tuples::element<N,T>::type type;
  };
}
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
