//----------------------------  extract_boundary_mesh_03.cc  ---------------------------
//    $Id: bem.cc 22693 2010-11-11 20:11:47Z kanschat $
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  extract_boundary_mesh_03.cc  ---------------------------


#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

using namespace std;


void test ()
{
  const unsigned int dim=2;

  Triangulation<dim-1,dim> boundary_mesh;
  map<Triangulation<dim-1,dim>::cell_iterator,Triangulation<dim,dim>::face_iterator >
    surface_to_volume_mapping;
  Triangulation<dim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);

  GridTools::extract_boundary_mesh (volume_mesh, boundary_mesh,
				    surface_to_volume_mapping);

  FE_Q <dim-1,dim>  boundary_fe (1);
  DoFHandler<dim-1,dim> boundary_dh(boundary_mesh);
  boundary_dh.distribute_dofs (boundary_fe);

  deallog << "n_dofs=" << boundary_dh.n_dofs() << std::endl;

  for (DoFHandler<dim-1,dim>::active_cell_iterator
	 cell = boundary_dh.begin_active(),
	 endc = boundary_dh.end(); cell!=endc; ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
	{
	  unsigned int index = cell->vertex_dof_index(v,0);
	  deallog << "vertex: " << v
		  << ", global: " << cell->vertex_index(v)
		  << " index: " << index << std::endl;
	}

      deallog << std::endl;
    }
}



int main ()
{

  ofstream logfile("extract_boundary_mesh_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test ();
 }
