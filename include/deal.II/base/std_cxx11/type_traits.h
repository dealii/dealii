// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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

#ifndef __deal2__std_cxx11_type_traits_h
#define __deal2__std_cxx11_type_traits_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#  include <type_traits>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
                                // TODO: could fill up with more types from
                                // C++11 type traits
  using std::is_pod;
  using std::is_standard_layout;
  using std::is_trivial;
}
DEAL_II_NAMESPACE_CLOSE

#else

#include <boost/type_traits.hpp>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using boost::is_pod;

                                // boost does not have is_standard_layout and
                                // is_trivial, but those are both a subset of
                                // is_pod
  template <typename T>
  struct is_standard_layout
  {
    static const bool value = boost::is_pod<T>::value;
  };

  template <typename T>
  struct is_trivial
  {
    static const bool value = boost::has_trivial_copy<T>::value &&
      boost::has_trivial_assign<T>::value &&
      boost::has_trivial_constructor<T>::value &&
      boost::has_trivial_destructor<T>::value;
  };
}
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
