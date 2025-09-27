// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
