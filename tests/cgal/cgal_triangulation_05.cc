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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/io.h>
#include <deal.II/cgal/triangulation.h>

#include "../tests.h"

using namespace CGALWrappers;

using K                 = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALTriangulation = CGAL::Delaunay_triangulation_3<K>;

template <int dim, int spacedim>
void
test()
{
  // Build a deal.II triangulation
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  // Use its vertices to build a Delaunay_triangulation_3
  CGALTriangulation tr;
  add_points_to_cgal_triangulation(tria.get_vertices(), tr);

  // Transform the Dealaunay_triangulation_3 to a deal.II triangulation
  tria.clear();
  cgal_triangulation_to_dealii_triangulation(tr, tria);
  deallog << "dim " << dim << ", spacedim " << spacedim << std::endl;
  GridOut go;
  go.write_vtk(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();
  test<3, 3>();
}
