//----------------------------  anisotropic_crash.cc  ---------------------------
//    anisotropic_crash.cc,v 1.1 2003/06/09 15:59:07 wolf Exp
//    Version:
//
//    Copyright (C) 2012 by the deal.II authors and Abner Salgado
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  anisotropic_crash.cc  ---------------------------


// Trying to catch a bug in the construction of patches when
// there is anisotropic refinement


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
  std::ofstream logfile ("anisotropic_crash/output");
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
  for( ; cell != end; ++cell )
    {
      switch (rand()%4)
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
  for( unsigned v=0; v<tri.n_used_vertices(); ++v )
    if (tri.get_used_vertices()[v] == true)
      GridTools::find_cells_adjacent_to_vertex( tri, v );
}
