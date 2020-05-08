// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// testing parent/ancestor and child relationship of CellIds

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

  deallog << std::boolalpha;

  for (const auto &pair :
       std::vector<std::pair<CellId, CellId>>{
         {id0, id0}, // same cell (expected: false, false)
         {id0, id1}, // same level (false, false)
         {id0, id2}, // child (true, true)
         {id0, id3}, // grand child (false, true)
         {id0, id4}} // cell with different coarse-cell id (false, false)
  )
    deallog << pair.first.is_parent_of(pair.second) << " "
            << pair.first.is_ancestor_of(pair.second) << std::endl;
}
