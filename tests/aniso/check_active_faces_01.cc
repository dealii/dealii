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

// n_faces() gives wrong numbers if anisotropic refinement is in place

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
testActiveFaceIndex(Triangulation<3> &tria)
{
  deallog << tria.n_active_cells() << " " << tria.n_cells() << " "
          << tria.n_hexs() << " " << tria.has_hanging_nodes() << std::endl;
  deallog << tria.n_active_faces() << " " << tria.n_faces() << " "
          << tria.n_quads() << " " << tria.has_hanging_nodes() << std::endl;
  std::vector<std::size_t> globalToActiveIndex(tria.n_quads(), 0);
  {
    std::size_t index = 0;
    for (auto &face : tria.active_face_iterators())
      {
        assert(!face->has_children());
        deallog << "Active Face Index " << index++ << " face index "
                << face->index() << " / " << tria.n_faces() << std::endl;
      }
  }
}

void
cubeTest()
{
  Triangulation<3> tria, triaInit;
  GridGenerator::hyper_cube(triaInit);
  triaInit.begin_active()->set_refine_flag(RefinementCase<3>::cut_x);
  triaInit.execute_coarsening_and_refinement();
  GridGenerator::flatten_triangulation(triaInit, tria);
  for (auto cell : tria.active_cell_iterators())
    {
      if (cell->center()(0) < 0.5 - 1e-6)
        {
          cell->set_refine_flag(RefinementCase<3>::cut_yz);
        }
      else
        {
          cell->set_refine_flag(RefinementCase<3>::cut_z);
        }
    }
  tria.execute_coarsening_and_refinement();
  for (auto cell : tria.active_cell_iterators())
    {
      if (cell->center()(0) < 0.5 - 1e-6)
        {
          cell->set_refine_flag(RefinementCase<3>::cut_z);
        }
      else
        {
          cell->set_refine_flag(RefinementCase<3>::cut_yz);
        }
    }
  tria.execute_coarsening_and_refinement();
  testActiveFaceIndex(tria);
}

int
main()
{
  initlog();
  cubeTest();
}
