//----------------------------  normal_vectors_01.cc  ---------------------------
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
//----------------------------  normal_vectors_01.cc  ---------------------------

/*
  Asking for normal vectors in the codim-1 case led to aborts.
*/

#include "../tests.h"

#include <base/quadrature_lib.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>

#include <fe/fe_q.h>
#include <fe/mapping_q.h>
#include <fe/fe_values.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

template <int dim>
void test ()
{
  Triangulation<dim-1,dim> tria;
  std::map<typename Triangulation<dim-1,dim>::cell_iterator,
    typename Triangulation<dim,dim>::face_iterator> surface_to_volume_mapping;

  HyperBallBoundary<dim> boundary_description;
  Triangulation<dim> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  
  volume_mesh.set_boundary (1, boundary_description);
  volume_mesh.set_boundary (0, boundary_description);
  volume_mesh.refine_global (1);
  
  static HyperBallBoundary<dim-1,dim> surface_description;
  tria.set_boundary (1, surface_description);
  tria.set_boundary (0, surface_description);
  
  std::set<unsigned char> boundary_ids;
  boundary_ids.insert(0);
  
  GridTools::extract_boundary_mesh (volume_mesh, tria,
				    surface_to_volume_mapping,
				    boundary_ids);
  
  FE_Q<dim-1,dim> fe (1);
  DoFHandler<dim-1,dim> dh(tria);
  MappingQ<dim-1,dim> mapping(1,true);
  
  dh.distribute_dofs (fe);
  
  FEFaceValues<dim-1,dim> fe_face_values (mapping, fe, QGauss<dim-2>(2), 
					  update_values         | update_quadrature_points  |
					  update_normal_vectors | update_JxW_values);
  
  typename DoFHandler<dim-1,dim>::active_cell_iterator
    cell = dh.begin_active();
  
  fe_face_values.reinit (cell,0);
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
