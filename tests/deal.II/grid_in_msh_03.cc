//----------------------------  grid_in_msh_03.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in_msh_03.cc  ---------------------------

// read a file in the MSH format used by the GMSH program. this
// particular mesh has type-15 cells (nodes) with more than one
// associated vertex. we failed to read the vertices

#include "../tests.h"
#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/grid_generator.h>
#include <base/logstream.h>

#include <fstream>
#include <iomanip>
#include <string>

std::ofstream logfile("grid_in_msh_03/output");


template<int dim>
void check_file (const std::string name,
		 typename GridIn<dim>::Format format)
{
  Triangulation<dim> tria;
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  gi.read(name, format);
  deallog << name
	  << '\t' << tria.n_vertices()
	  << '\t' << tria.n_cells()
	  << std::endl;

  GridOut grid_out;
  grid_out.write_gnuplot (tria, deallog.get_file_stream());
}

void filename_resolution()
{
  check_file<2> (std::string("grid_in_msh_03/mesh"), GridIn<2>::msh);
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  filename_resolution();
}

