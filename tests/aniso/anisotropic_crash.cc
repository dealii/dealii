// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// GridTools::find_cells_adjacent_to_vertex had a problem in that it
// wasn't prepared to deal with anisotropic refinement

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


/// to generate random numbers


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);
  deallog.get_file_stream().setf(std::ios::fixed);

  // Create triangulation
  Triangulation<2> tri;
  GridGenerator::hyper_cube(tri);
  tri.refine_global(2);

  // now do some sort of random, anisotropic
  // refinement
  Triangulation<2>::active_cell_iterator cell = tri.begin_active(),
                                         end  = tri.end();
  for (; cell != end; ++cell)
    {
      switch (Testing::rand() % 4)
        {
          /// If a randomly drawn
          /// number is 0 or 1 we
          /// cut x or y, resp.
          case 0:
            cell->set_refine_flag(RefinementCase<2>::cut_axis(0));
            break;
          case 1:
            cell->set_refine_flag(RefinementCase<2>::cut_axis(1));
            break;

          case 2:
            /// If the number is 2 we
            /// refine isotropically
            cell->set_refine_flag();
            break;

          default:
            /// If the number is 3 we don't refine
            ;
        }
    }
  /// refine the mesh
  tri.execute_coarsening_and_refinement();

  /// For each vertex find the patch of cells
  /// that surrounds it
  for (unsigned v = 0; v < tri.n_vertices(); ++v)
    if (tri.get_used_vertices()[v] == true)
      {
        deallog << "Vertex=" << v << std::endl;

        const std::vector<Triangulation<2>::active_cell_iterator> tmp =
          GridTools::find_cells_adjacent_to_vertex(tri, v);

        for (unsigned int i = 0; i < tmp.size(); ++i)
          deallog << "   " << tmp[i] << std::endl;
      }
}
