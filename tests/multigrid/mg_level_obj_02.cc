// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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
