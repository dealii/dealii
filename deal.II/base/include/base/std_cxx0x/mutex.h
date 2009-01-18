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
#ifndef __deal2__std_cxx0x_mutex_h
#define __deal2__std_cxx0x_mutex_h


#include <base/config.h>

#ifdef DEAL_II_CAN_USE_CXX0X

#  include <mutex>

#else

#  include <boost/thread/mutex.hpp>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx0x
{
  using boost::mutex;
}
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
