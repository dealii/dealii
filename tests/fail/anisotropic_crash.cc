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

/// to generate random numbers
#include <cstdlib>

using namespace dealii;

int main(){
				   /// Create triangulation, DOF handler and finite element space
  Triangulation<2> tri;
  DoFHandler<2> dh( tri );
  FE_Q<2> fe(1);
  GridGenerator::hyper_cube( tri );
  tri.refine_global( 2 );
  const unsigned n_cycles = 8;
  for( unsigned i=0; i<n_cycles; ++i ){
				     /// Construct the finite element space
    dh.distribute_dofs( fe );
				     /// For each vertex find the patch
    const unsigned n_dofs = dh.n_dofs(), n_vert = tri.n_used_vertices();
    std::vector< Point<2> > vertices = tri.get_vertices();
    for( unsigned v=0; v< n_vert; ++v ){
      std::vector< DoFHandler<2>::active_cell_iterator > patch
	= GridTools::find_cells_adjacent_to_vertex( dh, v );
    }
				     /// loop over cells and randomly refine
    DoFHandler<2>::active_cell_iterator cell = dh.begin_active(), end = dh.end();
    for( ; cell != end; ++cell ){
      int random = rand()%4;
				       /// If random is 0 or 1 we cut x or y, resp.
      if( random < 2 )
	cell->set_refine_flag( RefinementCase<2>::cut_axis( random ) );
				       /// If random is 2 we refine isotropically
      if( random == 2 )
	cell->set_refine_flag();
				       /// If random is 3 we don't refine
    }
				     /// refine the mesh
    tri.prepare_coarsening_and_refinement();
    tri.execute_coarsening_and_refinement();
  }
  dh.clear();
  return 0;
}
