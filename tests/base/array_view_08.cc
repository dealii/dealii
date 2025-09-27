// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test ArrayView::begin and ArrayView::end

#include <deal.II/base/array_view.h>

#include <type_traits>

#include "../tests.h"

template <typename T>
constexpr bool
is_const_reference()
{
  return std::is_reference_v<T> && std::is_const_v<std::remove_reference_t<T>>;
}

void
test()
{
  // test for a mutable view of a mutable vector
  {
    std::vector<int> v(10);
    ArrayView<int>   a(&v[4], 3);
    AssertThrow(a.begin() == &v[4], ExcInternalError());
    AssertThrow(a.end() == &v[7], ExcInternalError());

    static_assert(std::is_reference_v<decltype(*a.begin())>,
                  "type should be a reference");
    static_assert(std::is_reference_v<decltype(*a.end())>,
                  "type should be a reference");
    static_assert(!is_const_reference<decltype(*a.begin())>(),
                  "type should not be const");
    static_assert(!is_const_reference<decltype(*a.end())>(),
                  "type should not be const");
  }

  // and an immutable view of a mutable vector
  {
    std::vector<int>           v(10);
    const ArrayView<const int> a(&v[4], 3);
    AssertThrow(a.begin() == &v[4], ExcInternalError());
    AssertThrow(a.end() == &v[7], ExcInternalError());

    static_assert(is_const_reference<decltype(*a.begin())>(),
                  "type should be const reference");
    static_assert(is_const_reference<decltype(*a.end())>(),
                  "type should be const reference");
  }

  // and an immutable view of an immutable vector
  {
    const std::vector<int>     v(10, 42);
    const ArrayView<const int> a(&v[4], 3);
    AssertThrow(a.begin() == &v[4], ExcInternalError());
    AssertThrow(a.end() == &v[7], ExcInternalError());

    static_assert(is_const_reference<decltype(*a.begin())>(),
                  "type should be const reference");
    static_assert(is_const_reference<decltype(*a.end())>(),
                  "type should be const reference");
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
