// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
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


// verify that we can do things like cell->face() in 1d as well. here:
// Triangulation::get_boundary_ids() should return an empty vector when
// called to create a mesh that is a loop


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
test()
{
  Triangulation<2, 2> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);

  Triangulation<1, 2> tria;

  GridGenerator::extract_boundary_mesh(volume_mesh, tria);

  deallog << "n_cells = " << tria.n_active_cells() << std::endl;
  deallog << "n_boundary_ids = " << tria.get_boundary_ids().size() << std::endl;
}



int
main()
{
  initlog();

  test();

  return 0;
}
