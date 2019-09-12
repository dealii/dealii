// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// testing parent and child relationship of CellIds

#include <deal.II/grid/cell_id.h>

#include "../tests.h"

int
main()
{
  initlog();

  // create some arbitrary CellIds
  CellId id0(0, {0});
  CellId id1(0, {1});
  CellId id2(0, {0, 0});
  CellId id3(0, {0, 0, 0});
  CellId id4(1, {0});

  // same cell (expected: false)
  deallog << std::boolalpha;
  deallog << id0.is_parent_of(id0) << std::endl;

  // same level (false)
  deallog << id0.is_parent_of(id1) << std::endl;

  // child (true)
  deallog << id0.is_parent_of(id2) << std::endl;

  // grand child (false)
  deallog << id0.is_parent_of(id3) << std::endl;

  // cell with different coarse-cell id (false)
  deallog << id0.is_parent_of(id4) << std::endl;
}
