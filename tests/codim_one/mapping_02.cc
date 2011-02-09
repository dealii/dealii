//----------------------------  mapping_02.cc  ---------------------------
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
//----------------------------  mapping_02.cc  ---------------------------


// like _01, but use a quadratic mapping. since we now map line segments to
// curves, the normal vectors at different quadrature points should no longer
// be parallel

#include "../tests.h"

#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>
#include <fe/fe_values.h>
#include <fe/fe_q.h>
#include <fe/mapping_q.h>


template <int dim>
void test ()
{
  deallog << "Testing hyper_ball in dim: " << dim << "..."<< std::endl;

  Triangulation<dim> volume_mesh;
  GridGenerator::hyper_ball(volume_mesh);

  const HyperBallBoundary<dim-1,dim> surface_description;
  Triangulation<dim-1,dim> boundary_mesh;
  boundary_mesh.set_boundary (0, surface_description);

  GridTools::extract_boundary_mesh (volume_mesh, boundary_mesh);

  QGauss<dim-1> quadrature(dim == 2 ? 3 : 2);
  MappingQ<dim-1,dim> mapping(2);
  FE_Q<dim-1,dim> fe (1);
  
  FEValues<dim-1,dim> fe_values (mapping, fe, quadrature, update_cell_normal_vectors);
  
  for (typename Triangulation<dim-1,dim>::active_cell_iterator
	 cell = boundary_mesh.begin_active(); cell != boundary_mesh.end();
       ++cell)
    {
      deallog << "Cell = " << cell
	      << ", with center at " << cell->center()
	      << std::endl;
      fe_values.reinit (cell);

      for (unsigned int q=0; q<quadrature.size(); ++q)
	deallog << "  cell_normal[" << q << "] = "
		<< fe_values.cell_normal_vector(q)
		<< std::endl;
    }
}

  

int main ()
{
  std::ofstream logfile("mapping_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> ();
  test<3> ();

  return 0;
}
