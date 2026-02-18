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


// check MGLevelObject::clear()

#include <deal.II/base/mg_level_object.h>

#include <algorithm>

#include "../tests.h"

int
main()
{
  initlog();

  MGLevelObject<double> o(2, 4);

  AssertDimension(o.min_level(), 2);
  AssertDimension(o.max_level(), 4);
  AssertDimension(o.n_levels(), 3);

  o.clear();

  AssertDimension(o.min_level(), 0);
  AssertDimension(o.n_levels(), 0);

  deallog << "OK!" << std::endl;
}
