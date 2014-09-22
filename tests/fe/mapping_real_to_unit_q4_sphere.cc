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



// this is a variation of the _q4_curved test, with data extracted by
// Eric Heien from the Aspect code. we used to trip an assertion (the
// Newton iteration in transform_real_to_unit_cell did not converge)
// because the point we look for is outside the Q1 mapped cell but
// inside the Q4 mapped cell. as part of the algorithm, we get an
// initial guess for the inverse Q4 mapping by running an inverse Q1
// mapping, and we failed to find a solution in the latter.
//
// todo: this test used to fail before r30450 but the new mapping support
// points seem to fix this
//
// a redux of this test is in the mapping_real_to_unit_q4_sphere_*
// tests

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/fe/mapping_q.h>


void test_real_to_unit_cell()
{
  const unsigned int dim = 3;

  // define a spherical cap boundary
  // to be used for one of the faces
  const double radius = Point<dim>(1.43757e-10, 4.48023e+06, 4.48023e+06).norm();
  HyperBallBoundary<dim> boundary (Point<dim>(), radius);

  // create the mesh: a single cell
  // with the following coordinates:
  std::vector<Point<dim> > vertices;
  vertices.push_back (Point<dim>( 6.70384e-11, 3.17728e+06, 3.17728e+06));
  vertices.push_back (Point<dim>( -1.46060e+06, 3.99043e+06, 1.46060e+06));
  vertices.push_back (Point<dim>( -1.46060e+06, 1.46060e+06, 3.99043e+06));
  vertices.push_back (Point<dim>( -2.59424e+06, 2.59424e+06, 2.59424e+06));
  vertices.push_back (Point<dim>( 1.43757e-10, 4.48023e+06, 4.48023e+06));
  vertices.push_back (Point<dim>( -2.05956e+06, 5.62684e+06, 2.05956e+06));
  vertices.push_back (Point<dim>( -2.05956e+06, 2.05956e+06, 5.62684e+06));
  vertices.push_back (Point<dim>( -3.65809e+06, 3.65809e+06, 3.65809e+06));
  // the points above don't show
  // enough digits to have the same
  // outer radius, so normalize the
  // four outer ones
  for (unsigned int v=4; v<8; ++v)
    vertices[v] *= radius/vertices[v].norm();
  std::vector<CellData<dim> > cells;
  {
    CellData<dim> d;
    for (unsigned int i=0; i<8; ++i)
      d.vertices[i] = i;
    cells.push_back(d);
  }
  Triangulation<dim> triangulation;
  triangulation.create_triangulation (vertices, cells,
                                      SubCellData());

  // set the boundary indicator for
  // one face and adjacent edges of
  // the single cell
  triangulation.set_boundary (1, boundary);
  triangulation.begin_active()->face(5)->set_all_boundary_indicators (1);

  // now try to find the coordinates
  // of the following point in the
  // reference coordinate system of
  // the cell
  const Point<dim> p (-3.56413e+06, 1.74215e+06, 2.14762e+06);
  MappingQ<dim> map(4);
  Triangulation<dim >::active_cell_iterator
  cell = triangulation.begin_active();

  try
    {
      map.transform_real_to_unit_cell(cell,p);
    }
  catch (const Mapping<dim>::ExcTransformationFailed &)
    {
      deallog << "Point is outside!" << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_real_to_unit_cell();

  return 0;
}



