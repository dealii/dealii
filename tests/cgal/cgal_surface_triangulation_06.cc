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

// Convert a CGAL::Surface_mesh or a CGAL::Polyhedron_3 to a deal.II
// triangulation. CGAL surfaces are generated using an .off file describing
// the shape of the number '8'.

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/io.h>
#include <deal.II/cgal/triangulation.h>

#include <fstream>

#include "../tests.h"

using namespace CGALWrappers;
using CGALPoint  = CGAL::Point_3<CGAL::Simple_cartesian<double>>;
using Mesh       = CGAL::Surface_mesh<CGALPoint>;
using Polyhedron = CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>>;

int
main()
{
  initlog();
  Mesh          sm;
  std::ifstream input_sm(SOURCE_DIR "/input_grids/eight.off");
  input_sm >> sm;

  Triangulation<2, 3> tria;
  cgal_surface_mesh_to_dealii_triangulation(sm, tria);
  GridOut go;
  // std::ofstream out("surface_tria.vtk");
  // go.write_vtk(tria, out);
  go.write_vtk(tria, deallog.get_file_stream());

  // Now do the same with a Polyhedron
  deallog << "Polyhedron test" << std::endl;
  Polyhedron    P;
  std::ifstream input_p(SOURCE_DIR "/input_grids/eight.off");
  input_p >> P;


  tria.clear();
  cgal_surface_mesh_to_dealii_triangulation(P, tria);
  // std::ofstream output("poly_tria.vtk");
  // go.write_vtk(tria, output);
  go.write_vtk(tria, deallog.get_file_stream());
}
