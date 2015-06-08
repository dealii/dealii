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

#ifndef dealii__std_cxx11_unique_ptr_h
#define dealii__std_cxx11_unique_ptr_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CXX11

#  include <memory>
DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  using std::unique_ptr;
}
DEAL_II_NAMESPACE_CLOSE

#else

#include <boost/scoped_ptr.hpp>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx11
{
  /**
   * Implementation of a basic replacement for C++11's std::unique_ptr class.
   *
   * BOOST does not have a replacement for std::unique_ptr (because unique_ptr
   * requires move semantics that aren't available unless you have a C++11
   * compiler -- in which case you also have std::unique_ptr; see for example
   * http://stackoverflow.com/questions/2953530/unique-ptr-boost-equivalent)
   *
   * Consequently, we emulate the class by just wrapping a boost::scoped_ptr
   * in the cheapest possible way -- by just deriving from it and repeating
   * the basic constructors. Everything else is inherited from the scoped_ptr
   * class.
   *
   * There is no overhead to this approach: scoped_ptr cannot be copied or
   * moved. Instances of unique_ptr cannot be copied, and if you do not have a
   * C++11 compiler, then you cannot move anything anyway.
   */
  template <typename T>
  class unique_ptr : public boost::scoped_ptr<T>
  {
  public:
    unique_ptr () {}

    template<class Y>
    explicit unique_ptr (Y *p)
      :
      boost::scoped_ptr<T>(p)
    {}
  };

}
DEAL_II_NAMESPACE_CLOSE

#endif

// then allow using the old namespace name instead of the new one
DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x = std_cxx11;
DEAL_II_NAMESPACE_CLOSE

#endif
