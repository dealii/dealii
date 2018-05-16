// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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


// generate and refine a hyper shell in 3d with 6 cells. The other tests of
// this series test the other possibilities for the number of cells

#include "../tests.h"
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <iostream>

std::ofstream
logfile("output");


template <int dim>
void
check (double r1, double r2, unsigned int n)
{
  Point<dim> center;
  Triangulation<dim> tria (Triangulation<dim>::none);
  GridGenerator::hyper_shell (tria, center, r1, r2, n);
  static const SphericalManifold<dim> boundary(center);
  tria.set_manifold(0, boundary);
  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_gnuplot (tria, logfile);
}


int
main()
{
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  check<3> (.5, 1, 6);
}
