//----------------------------  cone_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  cone_01.cc  ---------------------------


// check ConeBoundary and GridGenerator::truncated_cone

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_c1.h>

#include <fstream>



template <int dim>
void check ()
{
  Triangulation<dim> triangulation;
  GridGenerator::truncated_cone (triangulation);
  Point<dim> p1, p2;
  p1[0] = -1;
  p2[0] = 1;
  static const ConeBoundary<dim> boundary (1, 0.5, p1, p2);
  triangulation.set_boundary (0, boundary);

  triangulation.refine_global (2);

  GridOut().write_gnuplot (triangulation,
			   deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("cone_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}



