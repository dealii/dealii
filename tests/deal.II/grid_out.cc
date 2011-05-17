//----------------------------  constraints.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  constraints.cc  ---------------------------


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


std::ofstream logfile("grid_out/output");


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

