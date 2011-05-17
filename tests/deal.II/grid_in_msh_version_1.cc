//----------------------------  grid_in_msh_version_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in_msh_version_1.cc  ---------------------------

// in implementing GMSH version 2 input, we forgot reading vertices in
// version 2 format as well. A fix by Victor Prosolin helped this
// problem

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>
#include <string>

std::ofstream logfile("grid_in_msh_version_1/output");


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
  check_file<2> (std::string("grid_in_msh_version_1/input_v1"), GridIn<2>::msh);
  check_file<2> (std::string("grid_in_msh_version_1/input_v2"), GridIn<2>::msh);
}


int main ()
{
  deallog << std::setprecision (5);
  logfile << std::setprecision (5);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  filename_resolution();
}

