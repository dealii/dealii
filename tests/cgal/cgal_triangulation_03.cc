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

// Convert a cgal Triangulation_2 to a dealii::Triangulation<2, ...>

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/io.h>
#include <CGAL/Triangulation_2.h>
#include <deal.II/cgal/triangulation.h>

#include "../tests.h"

using namespace CGALWrappers;

using K                 = CGAL::Simple_cartesian<double>;
using CGALTriangulation = CGAL::Triangulation_2<K>;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  CGALTriangulation tr;
  add_points_to_cgal_triangulation(tria.get_vertices(), tr);

  tria.clear();
  // Test the other way around
  cgal_triangulation_to_dealii_triangulation(tr, tria);
  deallog << "dim " << dim << ", spacedim " << spacedim << std::endl;
  GridOut go;
  go.write_vtk(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
}
