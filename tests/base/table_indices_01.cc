// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check TableIndices in various ways

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  const TableIndices<2> t1(84, 42);
  TableIndices<2>       t2;
  t2[0] = 84;
  t2[1] = 42;

  AssertThrow(t1 == t2, ExcInternalError());
  AssertThrow(t1[0] == t2[0], ExcInternalError());
  AssertThrow(t1[1] == t2[1], ExcInternalError());

  AssertThrow(!(t1 != t2), ExcInternalError());

  t2.sort();
  AssertThrow(t1 != t2, ExcInternalError());
  AssertThrow(t1[0] == t2[1], ExcInternalError());
  AssertThrow(t1[1] == t2[0], ExcInternalError());

  deallog << "OK" << std::endl;
}
