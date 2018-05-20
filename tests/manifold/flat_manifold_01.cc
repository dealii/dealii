// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

// Test that the flat manifold does what it should

#include "../tests.h"

// all include files you need here
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  deallog << "Testing dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  typename Triangulation<dim, spacedim>::active_cell_iterator cell;

  for(cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      // check that FlatManifold returns the middle of the cell.
      deallog << "Cell: " << cell << std::endl;
      if(cell->get_manifold().get_new_point_on_cell(cell).distance(
           cell->center())
         > 1e-6)
        {
          deallog << "Default manifold: "
                  << cell->get_manifold().get_new_point_on_cell(cell)
                  << std::endl;
          deallog << "Center of cell  : " << cell->center() << std::endl;
        }
      else
        {
          deallog << "OK!" << std::endl;
        }
    }
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
