// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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


// Create a hyper ball, refine it, extract an iges of the boundary,
// and create a new mesh using gmsh from that iges file.

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/gmsh/utilities.h>
#include <deal.II/opencascade/utilities.h>
#include <deal.II/opencascade/boundary_lib.h>

int main ()
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::hyper_ball(tria);

  tria.refine_global(4);

  auto curves = OpenCASCADE::create_curves_from_triangulation_boundary(tria);
  tria.clear();

  Gmsh::create_triangulation_from_boundary_curve(curves[0], tria);

  GridOut go;
  std::ofstream ofile("output.svg");
  go.write_svg(tria, ofile);
  go.write_svg(tria, deallog.get_file_stream());

  return 0;
}
