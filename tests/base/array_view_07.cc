// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test make_const_array_view

#include <deal.II/base/array_view.h>

#include "../tests.h"

template <typename ElementType>
void
const_foo(const ArrayView<const std::vector<ElementType>> &view)
{
  AssertThrow(view[0][0] == 1, ExcInternalError());
  deallog << "OK" << std::endl;
}

template <typename ElementType>
void
foo(const ArrayView<std::vector<ElementType>> &view)
{
  AssertThrow(view[0][0] == 1, ExcInternalError());
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  std::vector<std::vector<int>> v(1, std::vector<int>(1, 1));
  // this doesn't work
  // const_foo(make_array_view(v));
  const_foo(make_const_array_view(v));
  foo(make_array_view(v));
  // this doesn't work
  // foo(make_const_array_view(v));
}
