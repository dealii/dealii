// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

// Test that the flat manifold does what it should. This time on faces.

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

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      // check that FlatManifold returns the middle of the cell.
      deallog << "Cell: " << cell << std::endl;
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        {
          const typename Triangulation<dim, spacedim>::face_iterator &face =
            cell->face(f);
          if (face->get_manifold().get_new_point_on_face(face).distance(
                face->center()) > 1e-6)
            {
              deallog << "Face            : " << face << std::endl;
              deallog << "Default manifold: "
                      << cell->get_manifold().get_new_point_on_cell(cell)
                      << std::endl;
              deallog << "Center of cell  : " << cell->center() << std::endl;
            }
          else
            {
              deallog << "Face " << face << " is OK!" << std::endl;
            }
        }
    }
}

int
main()
{
  initlog();

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
