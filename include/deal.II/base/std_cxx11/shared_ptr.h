// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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

#ifndef dealii__std_cxx11_shared_ptr_h
#define dealii__std_cxx11_shared_ptr_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#  include <memory>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using std::shared_ptr;
  using std::enable_shared_from_this;
  using std::make_shared;
  using std::dynamic_pointer_cast;
}
DEAL_II_NAMESPACE_CLOSE

#else

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/make_shared.hpp>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using boost::shared_ptr;
  using boost::enable_shared_from_this;
  using boost::make_shared;
  using boost::dynamic_pointer_cast;
}
DEAL_II_NAMESPACE_CLOSE

#endif

// then allow using the old namespace name instead of the new one
DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x = std_cxx11;
DEAL_II_NAMESPACE_CLOSE

#endif
