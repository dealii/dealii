// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

// take a 3d mesh, take all vertices, shift them a little bit and check that
// we correctly identify the closest vertex position
// The result should be an increasing sequence of numbers

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

void check(Triangulation<3>& tria)
{
  const std::vector<Point<3>>& v = tria.get_vertices();
  for(unsigned i = 0; i < v.size(); i++)
    deallog << "["
            << GridTools::find_closest_vertex(
                 tria, v[i] + Point<3>(0.01, -0.01, 0.01))
            << "] ";

  deallog << std::endl;
}

int
main()
{
  initlog();

  try
    {
      Triangulation<3> coarse_grid;
      GridGenerator::hyper_cube(coarse_grid);
      coarse_grid.refine_global(3);
      check(coarse_grid);
    }
  catch(const std::exception& exc)
    {
      // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}
