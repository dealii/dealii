// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// testing serialize function for class CellId

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int spacedim>
void
test(const unsigned int &ref)
{
  Triangulation<spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(ref);
  std::vector<CellId> cell_ids;
  for (auto &cell : tria.active_cell_iterators())
    {
      cell_ids.push_back(cell->id());
    }

  auto buffer = Utilities::pack(cell_ids);

  auto unpacked = Utilities::unpack<std::vector<CellId>>(buffer);

  unsigned int i  = 0;
  bool         ok = true;
  for (auto &cell : tria.active_cell_iterators())
    {
      if (cell->id() != unpacked[i++])
        {
          deallog << "NOT OK; problem with cell " << i << std::endl;
          ok = false;
        }
    }

  if (ok)
    deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test<1>(4);
  test<2>(3);
  test<3>(2);
}
