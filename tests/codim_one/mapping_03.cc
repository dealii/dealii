// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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
  for (unsigned int i=0; i<dim; ++i)
    p[i] = .2;

  Point<dim-1> q =
    mapping.transform_real_to_unit_cell(mesh.begin_active(), p);

  deallog << "Mapping Q("<< degree<< "): P: " << p
          << ", on unit: " << q << std::endl;

}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> (1);
  test<2> (2);

  test<3> (1);
  test<3> (2);

  return 0;
}
