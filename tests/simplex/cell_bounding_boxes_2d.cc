// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Test bounding box creation for simplex 2D tria.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
create_simplex_grid(Triangulation<dim> &tria_simplex)
{
  tria_simplex.clear();
  Triangulation<dim> tria_standard;
  GridGenerator::hyper_cube(tria_standard);
  GridGenerator::convert_hypercube_to_simplex_mesh(tria_standard, tria_simplex);
}


template <int dim>
void
compute_bounding_boxes(Triangulation<dim> &tria)
{
  for (const auto cell : tria.active_cell_iterators())
    {
      const auto cell_bb = cell->bounding_box();
      deallog << "Center : " << cell_bb.center()
              << " , volume : " << cell_bb.volume() << std::endl;
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(5);
  Triangulation<2> tria;
  create_simplex_grid(tria);
  compute_bounding_boxes(tria);
}
