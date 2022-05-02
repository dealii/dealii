// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------
#ifndef dealii_cxx20_iota_view_h
#define dealii_cxx20_iota_view_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_HAVE_CXX20
#  include <ranges>
#else
#  include <boost/range/irange.hpp>
#endif

DEAL_II_NAMESPACE_OPEN

namespace std_cxx20
{
  namespace ranges
  {
#ifndef DEAL_II_HAVE_CXX20
    /**
     * A poor-man's implementation of std::ranges::iota_view using
     * boost's integer_range class. The two classes are not completely
     * identical, but close enough so that we can use the boost
     * class until std::ranges::iota_range becomes available in
     * compilers.
     *
     * The std::ranges::iota_view class represents a range of values
     * or iterators that one can loop over. A documentation of this
     * class can be found at
     * https://en.cppreference.com/w/cpp/ranges/iota_view .
     */
    template <typename IncrementableType, typename /*BoundType*/>
    using iota_view = boost::integer_range<IncrementableType>;
#else
    using std::ranges::iota_view;
#endif
  } // namespace ranges
} // namespace std_cxx20

DEAL_II_NAMESPACE_CLOSE

#endif
