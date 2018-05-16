// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
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


// GridTools::regularize_corner_cells on more complicated mesh

#include "../tests.h"
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

int
main()
{
  initlog();

  Point<2> p0(0,0), p1(4,1);
  Point<2> c0(.1, .5), c1(3.9,.5);

  SphericalManifold<2> m0(c0);
  SphericalManifold<2> m1(c1);

  Triangulation<2> tria;
  std::vector<unsigned int> subdivisions(2);
  subdivisions[0] = 4;
  subdivisions[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(tria,subdivisions, p0,p1,true);

  GridTools::copy_boundary_to_manifold_id(tria);

  tria.set_manifold(0, m0);
  tria.set_manifold(1, m1);

  GridTools::regularize_corner_cells (tria);
  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_msh (tria, deallog.get_file_stream());

  return 0;
}


