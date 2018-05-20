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

#ifndef dealii_std_cxx11_bind_h
#define dealii_std_cxx11_bind_h

#include <deal.II/base/config.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN
// In boost, the placeholders _1, _2, ... are in the global namespace. In
// C++11, they are in namespace std::placeholders, which makes them awkward to
// use. Import them into the dealii::std_cxx11 namespace instead and do the
// same below if we use boost instead. Namespace 'placeholders' is also defined
// in dealii::std_cxx11 namespace to make code C++ standard compatible.
// That is to say, if std::something works with C++11 standard,
// then dealii::std_cxx11::something should also work.
namespace std_cxx11
{
  using namespace std::placeholders;
  using std::bind;
  using std::cref;
  using std::ref;
  using std::reference_wrapper;

  namespace placeholders = std::placeholders;
} // namespace std_cxx11
DEAL_II_NAMESPACE_CLOSE

// then allow using the old namespace name instead of the new one
DEAL_II_NAMESPACE_OPEN
namespace std_cxx1x = std_cxx11;
DEAL_II_NAMESPACE_CLOSE

#endif
