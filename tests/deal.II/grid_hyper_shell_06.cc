// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
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


// GridGenerator::hyper_shell colorized the faces but forgot the
// edges. This is not useful because the colorization is usually done
// so that one can attach a boundary or manifold object to these parts
// of the boundary

#include "../tests.h"
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

std::ofstream logfile("output");


template<int dim>
void check (double r1, double r2, unsigned int n)
{
  Point<dim> center;
  Triangulation<dim> tria (Triangulation<dim>::none);
  GridGenerator::hyper_shell (tria, center, r1, r2, n, true);
  static const HyperShellBoundary<dim> boundary(center);
  tria.set_boundary(0, boundary);

  for (typename Triangulation<dim>::cell_iterator cell=tria.begin();
       cell != tria.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
	for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
	  Assert (cell->face(f)->line(l)->boundary_indicator()
		  ==
		  cell->face(f)->boundary_indicator(),
		  ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<3> (.5, 1, 6);
  check<3> (.5, 1, 12);
  check<3> (.5, 1, 96);
}
