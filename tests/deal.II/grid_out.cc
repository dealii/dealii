// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>


std::ofstream logfile("output");


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  static const HyperBallBoundary<dim> x;
  if (dim == 2)
    {
      tria.set_boundary (0, x);
      GridGenerator::hyper_ball (tria);
    }
  else
    GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  GridOut grid_out;
  GridOutFlags::Eps<2> eps2(GridOutFlags::EpsFlagsBase::width,
                            300, .5, false, 5, true);
  grid_out.set_flags (eps2);

  if (dim != 1)
    grid_out.write_eps (tria, logfile);
  grid_out.write_gnuplot (tria, logfile);
  grid_out.set_flags (GridOutFlags::Ucd(true));
  grid_out.write_ucd (tria, logfile);
  if (dim != 1)
    grid_out.write_dx (tria, logfile);
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
}

