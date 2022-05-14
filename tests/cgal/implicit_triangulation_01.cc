// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Convert a cgal Delaunay_triangulation_3 to a
// dealii::Triangulation<dim, spacedim>
// Non trivial case of a sphere. We need Exact predicates, inexact construction
// kernel here, since Simple_cartesian will throw an exception when trying to
// add some almost coplanar points.

#include <deal.II/base/config.h>

#include <deal.II/base/function_parser.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace CGAL::parameters;

int
main()
{
  initlog();
  // Build a deal.II triangulation
  Triangulation<3>  tria;
  FunctionParser<3> implicit_function("(1-sqrt(x^2+y^2))^2+z^2-.25");
  GridGenerator::implicit_function(
    tria, implicit_function, Point<3>(1, 0, 0), 10.0, cell_size = 0.4);
  GridOut go;
  {
    std::ofstream of("tria.vtk");
    go.write_vtk(tria, of);
  }
  remove("tria.vtk");
  // If we got here, everything was ok, including writing the grid.
  deallog << "OK" << std::endl;
}
