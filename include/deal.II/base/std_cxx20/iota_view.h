// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
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
