//----------------------------  constraints.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  constraints.cc  ---------------------------


#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_generator.h>
#include <base/logstream.h>

#include <fstream>


std::ofstream logfile("grid_out.output");


template <int dim>
void test ()
{  
  Triangulation<dim> tria;
  if (dim == 2)
    {
      static const HyperBallBoundary<dim> x;
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
  grid_out.write_ucd (tria, logfile);
  if (dim != 1)
    grid_out.write_dx (tria, logfile);
}


int main ()
{
  logfile.precision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
}

