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

#ifndef __deal2__std_cxx11_bind_h
#define __deal2__std_cxx11_bind_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#  include <functional>

DEAL_II_NAMESPACE_OPEN
// in boost, the placeholders _1, _2, ... are in the global namespace. in
// C++0x, they are in namespace std::placeholders, which makes them awkward to
// use. import them into the deal.II::std_cxx11 namespace instead and do them
// same below if we use boost instead.
namespace std_cxx11
{
  using namespace std::placeholders;
  using std::bind;
  using std::ref;
  using std::cref;
  using std::reference_wrapper;
}
DEAL_II_NAMESPACE_CLOSE

#else

#include <boost/bind.hpp>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using boost::bind;
  using boost::ref;
  using boost::cref;
  using boost::reference_wrapper;
  
  // now also import the _1, _2 placeholders from the global namespace
  // into the current one as suggested above
  using ::_1;
  using ::_2;
  using ::_3;
  using ::_4;
  using ::_5;
  using ::_6;
  using ::_7;
  using ::_8;
  using ::_9;
}
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
