// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Check neighborhood information on anisotropic grids when doing single
// direction refinement (only cut_x, cut_y, or cut_z)

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
print_cell_info(const Triangulation<3> &tria)
{
  for (const auto &cell : tria.cell_iterators())
    {
      deallog << "Cell : " << cell->level() << " | " << cell->index()
              << std::endl;
      for (const auto i : cell->face_indices())
        if (!cell->face(i)->at_boundary())
          deallog << "Neighbour : " << cell->neighbor(i)->level() << " | "
                  << cell->neighbor(i)->index() << std::endl;
    }
}

void
test()
{
  Triangulation<3> tria;
  GridGenerator::hyper_cube(tria);
  tria.begin_active()->set_refine_flag(RefinementCase<3>::cut_x);
  tria.execute_coarsening_and_refinement();
  for (const auto &cell : tria.active_cell_iterators())
    {
      cell->set_refine_flag(RefinementCase<3>::cut_y);
    }
  tria.execute_coarsening_and_refinement();

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center()(0) > 0.5)
      cell->set_refine_flag(RefinementCase<3>::cut_z);
  tria.execute_coarsening_and_refinement();

  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->center()(0) < 0.5)
        cell->set_refine_flag(RefinementCase<3>::cut_z);
    }
  tria.execute_coarsening_and_refinement();

  for (const auto &cell : tria.active_cell_iterators())
    cell->set_refine_flag(RefinementCase<3>::cut_y);
  tria.execute_coarsening_and_refinement();

  print_cell_info(tria);
}

int
main()
{
  initlog();
  test();
}
