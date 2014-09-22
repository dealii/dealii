// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// document bug in GridGenerator::truncated_cone
// (no cells are returned if half_length < 0.5*radius for dim=3)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_c1.h>

#include <fstream>



template <int dim>
void check ()
{
  deallog << "dim=" << dim << std::endl;
  
  Triangulation<dim> triangulation;
  double r1 = 0.5, r2 = 1.0, halfl = 0.25;
  GridGenerator::truncated_cone (triangulation, r1, r2, halfl);
  Point<dim> p1, p2;
  p1[0] = -halfl;
  p2[0] = halfl;
  static const ConeBoundary<dim> boundary (r1, r2, p1, p2);
  triangulation.set_boundary (0, boundary);

  triangulation.refine_global (2);

  GridOut().write_gnuplot (triangulation,
                           deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}



