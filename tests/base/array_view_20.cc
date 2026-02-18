// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test for class ArrayView

#include <deal.II/base/array_view.h>

#include "../tests.h"


void
test()
{
  std::vector<int> v(10);
  ArrayView<int>   a(v);
  Assert(!a.empty(), ExcInternalError());

  std::vector<int> x;
  ArrayView<int>   b(x);
  Assert(b.empty(), ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
