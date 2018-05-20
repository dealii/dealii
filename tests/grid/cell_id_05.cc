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

// testing serialize function for class CellId

#include "../tests.h"

#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

template <int spacedim>
void
test(const unsigned int& ref)
{
  Triangulation<spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(ref);
  std::vector<CellId> cell_ids;
  for(auto cell : tria.active_cell_iterators())
    {
      cell_ids.push_back(cell->id());
    }

  auto buffer = Utilities::pack(cell_ids);

  auto unpacked = Utilities::unpack<std::vector<CellId>>(buffer);

  unsigned int i  = 0;
  bool         ok = true;
  for(auto cell : tria.active_cell_iterators())
    {
      if(cell->id() != unpacked[i++])
        {
          deallog << "NOT OK; problem with cell " << i << std::endl;
          ok = false;
        }
    }

  if(ok)
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
