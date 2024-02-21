// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test refinement of simplex mesh.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <set>

#include "../tests.h"

#include "./simplex_grids.h"

void
test(const unsigned int v)
{
  const unsigned int dim = 2;

  Triangulation<dim> tria;

  if (v < 4)
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);
  else
    GridGenerator::subdivided_hyper_cube_with_simplices_mix(tria, 2);


  if (v == 1 || v == 4)
    {
      tria.begin()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
  else if (v == 2 || v == 5)
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

  // verify that all faces have consistently set orientation flags.
  std::set<std::pair<unsigned int, unsigned int>> all_faces;
  for (const auto &cell : tria.active_cell_iterators())
    {
      const auto v0 = cell->vertex_index(0);
      const auto v1 = cell->vertex_index(1);
      const auto v2 = cell->vertex_index(2);

      const auto p0  = cell->face_orientation(0) ? std::make_pair(v0, v1) :
                                                   std::make_pair(v1, v0);
      const auto p1  = cell->face_orientation(1) ? std::make_pair(v1, v2) :
                                                   std::make_pair(v2, v1);
      const auto p2  = cell->face_orientation(2) ? std::make_pair(v2, v0) :
                                                   std::make_pair(v0, v2);
      const auto p0w = std::make_pair(p0.second, p0.first);
      const auto p1w = std::make_pair(p1.second, p1.first);
      const auto p2w = std::make_pair(p2.second, p2.first);

      all_faces.insert(p0);
      all_faces.insert(p1);
      all_faces.insert(p2);

      if (all_faces.count(p0w) == 1)
        deallog << "found inconsistent line (" << p0w.first << ", "
                << p0w.second << ")" << std::endl;
      if (all_faces.count(p1w) == 1)
        deallog << "found inconsistent line (" << p1w.first << ", "
                << p1w.second << ")" << std::endl;
      if (all_faces.count(p2w) == 1)
        deallog << "found inconsistent line (" << p2w.first << ", "
                << p2w.second << ")" << std::endl;
    }

  GridOut grid_out;
#if false
  std::ofstream out("mesh." + std::to_string(v) + ".vtk");
  grid_out.write_vtk(tria, out);
#else
  grid_out.write_vtk(tria, deallog.get_file_stream());
#endif
}

int
main()
{
  initlog();
  test(0); // no refinement
  test(1); // refinement of a single cell
  test(2); // global refinement
  test(3); // refine left bottom corner twice
  test(4); // mixed mesh: refinement of a single cell
  test(5); // mixed mesh: global refinement
}
