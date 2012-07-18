//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__std_cxx1x_type_traits_h
#define __deal2__std_cxx1x_type_traits_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_CAN_USE_CXX1X

#  include <type_traits>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x
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
namespace std_cxx1x
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
