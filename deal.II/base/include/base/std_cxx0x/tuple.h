//---------------------------------------------------------------------------
//    $Id: std_cxx0x_tuple.h 16094 2008-05-14 19:18:25Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__std_cxx0x_tuple_h
#define __deal2__std_cxx0x_tuple_h


#include <base/config.h>

#ifdef DEAL_II_CAN_USE_CXX0X

#  include <tuple>

#else

#include <boost/tuple/tuple.hpp>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx0x
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
