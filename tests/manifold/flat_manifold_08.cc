// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// Test that Manifold::get_normal_vector also works for a special case that
// used to crash with another algorithm.

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
  std::vector<Point<spacedim>> vertices(GeometryInfo<dim>::vertices_per_cell);
  std::vector<CellData<dim>>   cells(1);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  vertices[0] = Point<dim>(0, 0, 0);
  vertices[1] = Point<dim>(1, 0, 0);
  vertices[2] = Point<dim>(0.5, 0.4, 0);
  vertices[3] = Point<dim>(1.5, 0.4, 0);

  vertices[4] = Point<dim>(0, 0, 1);
  vertices[5] = Point<dim>(1, 0, 1);
  vertices[6] = Point<dim>(0.5, 0.4, 1);
  vertices[7] = Point<dim>(1.5, 0.4, 1);

  Triangulation<dim, spacedim> tria;
  tria.create_triangulation(vertices, cells, SubCellData());

  typename Triangulation<dim, spacedim>::active_cell_iterator cell =
    tria.begin_active();

  Point<dim> p1(0.5, 0, 0);
  deallog << "Normal vector of face 4: "
          << cell->get_manifold().normal_vector(cell->face(4), p1) << std::endl;
  deallog << "Center of face 4: " << cell->face(4)->center() << std::endl;
}

int
main()
{
  initlog();

  test<3, 3>();

  return 0;
}
