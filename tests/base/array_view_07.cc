// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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
