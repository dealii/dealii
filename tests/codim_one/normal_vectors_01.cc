//----------------------------  normal_vectors_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  normal_vectors_01.cc  ---------------------------

//  Asking for face normal vectors in the codim-1 case led to aborts.


#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

template <int dim>
void test ()
{
  Triangulation<dim,dim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  
  Triangulation<dim-1,dim> tria;
  GridTools::extract_boundary_mesh (volume_mesh, tria);
  
  FE_Q<dim-1,dim> fe (1);
  DoFHandler<dim-1,dim> dh(tria);
  MappingQ<dim-1,dim> mapping(1,true);
  
  dh.distribute_dofs (fe);
  
  FEFaceValues<dim-1,dim> fe_face_values (mapping, fe, QTrapez<dim-2>(),
					  update_normal_vectors);

  for (typename DoFHandler<dim-1,dim>::active_cell_iterator
	 cell = dh.begin_active();
       cell != dh.end(); ++cell)
    {
      deallog << "Face centered at " << cell->center()
	      << std::endl;

      for (unsigned int f=0; f<GeometryInfo<dim-1>::faces_per_cell; ++f)
	{
	  deallog << "  Edge centered at " << cell->face(f)->center()
		  << std::endl;
      
	  fe_face_values.reinit (cell,f);
	  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
	    deallog << "    normal_vector=" << fe_face_values.normal_vector(q)
		    << std::endl;
	}
    }
}



int main ()
{
  std::ofstream logfile("normal_vectors_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

//  test<2> ();
  test<3> ();
  
  return 0;
}
