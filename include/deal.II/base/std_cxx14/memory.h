// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
#ifndef dealii_cxx14_memory_h
#define dealii_cxx14_memory_h

#include <deal.II/base/config.h>

#include <memory>

#ifndef DEAL_II_WITH_CXX14
// needed for array type check
#  include <type_traits>
#endif

DEAL_II_NAMESPACE_OPEN
namespace std_cxx14
{
#ifdef DEAL_II_WITH_CXX14
  using std::make_unique;
#else
  template <typename T, typename... Args>
  inline
  std::unique_ptr<T>
  make_unique(Args &&... constructor_arguments)
  {
    static_assert(!std::is_array<T>::value,
                  "This function is not implemented for array types.");
    return std::unique_ptr<T>(new T(std::forward<Args>(constructor_arguments)...));
  }
#endif
}
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx14_memory_h
