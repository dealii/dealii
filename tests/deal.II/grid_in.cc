//----------------------------  constraints.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002 by the deal.II authors
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
#include <grid/grid_in.h>
#include <grid/grid_generator.h>
#include <base/logstream.h>

#include <fstream>


std::ofstream logfile("grid_in.output");


template <int dim>
void test ()
{  
  Triangulation<dim> tria;
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in ("grid_in.in");
  gi.read_ucd (in);
  
  static const HyperBallBoundary<dim> x;
  tria.set_boundary (0, x);
//  tria.refine_global(1);
  
  GridOut grid_out;
  grid_out.write_ucd (tria, logfile);
};


int main ()
{
  logfile.precision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2> ();
};

