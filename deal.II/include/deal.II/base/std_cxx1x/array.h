//---------------------------------------------------------------------------
//    $Id$
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
#ifndef __deal2__std_cxx1x_array_h
#define __deal2__std_cxx1x_array_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_CAN_USE_CXX1X

#  include <array>

#else

#include <boost/array.hpp>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x
{
  using boost::array;
}
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
