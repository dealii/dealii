
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// test VectorTools::interpolate_boundary_values for codim=1. like _01
// but for 1d triangulations

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vectors.h>

#include <string>

std::ofstream logfile("interpolate_boundary_values_03/output");

void test() {
  const int dim = 1;
  const int spacedim = 2;
  
  Triangulation<dim, spacedim> tria;
  Triangulation<spacedim> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  std::set<unsigned char> boundary_ids;
  boundary_ids.insert(0);
  
  GridTools::extract_boundary_mesh (volume_mesh, tria,boundary_ids);
    
  deallog << tria.n_active_cells() << " active cells" << std::endl;

  FE_Q<dim,spacedim> fe(1);
  DoFHandler<dim,spacedim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;
  
				   // test left and right boundary
				   // separately
  for (unsigned int boundary_id=0; boundary_id<2; ++boundary_id)
    {
      std::map<unsigned int, double> bv;
      VectorTools::interpolate_boundary_values (dof_handler,
						boundary_id,
						Functions::SquareFunction<spacedim>(),
						bv);
      deallog << bv.size() << " boundary degrees of freedom" << std::endl;
      
      for (std::map<unsigned int, double>::const_iterator i = bv.begin();
	   i != bv.end(); ++i)
	deallog << i->first << ' ' << i->second << std::endl;
      
      for (DoFHandler<dim,spacedim>::active_cell_iterator
      	     cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      	  if (cell->at_boundary(f) &&
      	      (cell->face(f)->boundary_indicator() == boundary_id))
      	    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
      	      for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
      		{
      		  Assert (bv.find(cell->face(f)->vertex_dof_index(v,i))
      			  != bv.end(),
      			  ExcInternalError());
      		  Assert (bv[cell->face(f)->vertex_dof_index(v,i)]
      			  ==
      			  Functions::SquareFunction<spacedim>()
      			  .value(cell->face(f)->vertex(v),i),
      			  ExcInternalError());
      		}
    }
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  test();
  
  return 0;
}

