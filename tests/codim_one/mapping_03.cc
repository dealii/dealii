//----------------------------  mapping_02.cc  ---------------------------
//    $Id: bem.cc 22693 2010-11-11 20:11:47Z kanschat $
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>


template <int dim>
void test (unsigned int degree)
{
  Triangulation<dim-1, dim> mesh;
  GridGenerator::hyper_cube(mesh);

  QGauss<dim-1> quadrature(dim == 2 ? 3 : 2);
  MappingQ<dim-1,dim> mapping(degree);
  Point<dim> p;

				   // Try to project a point on the
				   // surface
  for(unsigned int i=0; i<dim; ++i)
    p[i] = .2;

  Point<dim-1> q =
    mapping.transform_real_to_unit_cell(mesh.begin_active(), p);

  deallog << "Mapping Q("<< degree<< "): P: " << p
	  << ", on unit: " << q << std::endl;

}



int main ()
{
  std::ofstream logfile("mapping_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> (1);
  test<2> (2);

  test<3> (1);
  test<3> (2);

  return 0;
}
