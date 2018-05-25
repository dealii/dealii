// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// test that is_contiguous works correctly for pointers

#include <deal.II/base/array_view.h>

#include "../tests.h"

void
test()
{
  constexpr int *p = nullptr;
  constexpr int *q = nullptr;

  // the following code would fail if we tried to call the
  // non-constexpr version of the function, so this really must be the
  // pointer overload which we know always returns true
  constexpr bool b = internal::ArrayViewHelper::is_contiguous(p, q);

  Assert(b == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();

  test();
}
