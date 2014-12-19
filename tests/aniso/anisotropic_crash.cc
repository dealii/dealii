// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



// GridTools::find_cells_adjacent_to_vertex had a problem in that it
// wasn't prepared to deal with anisotropic refinement

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/base/logstream.h>

/// to generate random numbers
#include <cstdlib>
#include <fstream>

using namespace dealii;

int main()
{
  std::ofstream logfile ("output");
  logfile.precision (3);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // Create triangulation
  Triangulation<2> tri;
  GridGenerator::hyper_cube( tri );
  tri.refine_global( 2 );

  // now do some sort of random, anisotropic
  // refinement
  Triangulation<2>::active_cell_iterator cell = tri.begin_active(), end = tri.end();
  for ( ; cell != end; ++cell )
    {
      switch (Testing::rand()%4)
        {
          /// If a randomly drawn
          /// number is 0 or 1 we
          /// cut x or y, resp.
        case 0:
          cell->set_refine_flag( RefinementCase<2>::cut_axis(0) );
          break;
        case 1:
          cell->set_refine_flag( RefinementCase<2>::cut_axis(1) );
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
  for ( unsigned v=0; v<tri.n_vertices(); ++v )
    if (tri.get_used_vertices()[v] == true)
      {
        deallog << "Vertex=" << v << std::endl;

        const std::vector<Triangulation<2>::active_cell_iterator>
        tmp = GridTools::find_cells_adjacent_to_vertex( tri, v );

        for (unsigned int i=0; i<tmp.size(); ++i)
          deallog << "   " << tmp[i] << std::endl;
      }
}
