// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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



// Test the correct setup of neighbors during refinement of simplex mesh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "./simplex_grids.h"

void
test(const unsigned int v)
{
  const unsigned int dim = 2;

  Triangulation<dim> tria;

  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  if (v == 1)
    {
      tria.begin()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
  else if (v == 2)
    {
      tria.refine_global(1);
    }
  else if (v == 3)
    {
      tria.begin()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      tria.begin_active(1)->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << cell->level() << ':' << cell->index() << "     ";
      for (const auto f : cell->face_indices())
        if (cell->at_boundary(f))
          deallog << "-:- ";
        else
          deallog << cell->neighbor_level(f) << ':' << cell->neighbor_index(f)
                  << ' ';
      deallog << std::endl;
    }

  deallog << std::endl;
}

int
main()
{
  initlog();
  test(0); // no refinement
  test(1); // refinement of a single cell
  test(2); // global refinement
}
