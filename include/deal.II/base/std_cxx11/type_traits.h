// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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

#ifndef dealii__std_cxx11_type_traits_h
#define dealii__std_cxx11_type_traits_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#  include <type_traits>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  // TODO: could fill up with more types from
  // C++11 type traits
  using std::is_fundamental;
  using std::is_pod;
  using std::is_pointer;
  using std::is_standard_layout;
  using std::is_trivial;
  using std::enable_if;
}
DEAL_II_NAMESPACE_CLOSE

#else

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/type_traits.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION<105600
#include <boost/utility/enable_if.hpp>
#else
#include <boost/core/enable_if.hpp>
#endif
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using boost::is_fundamental;
  using boost::is_pod;
  using boost::is_pointer;

  // boost::enable_if_c, *not* boost::enable_if, is equivalent to std::enable_if.
  template <bool B, class T = void>
  struct enable_if : public boost::enable_if_c<B, T>
  {};

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


// then allow using the old namespace name instead of the new one
DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x = std_cxx11;
DEAL_II_NAMESPACE_CLOSE

#endif
