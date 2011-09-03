//----------------------------  grid_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_tools.cc  ---------------------------


// check GridTools::diameter for codim-1 meshes


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>

std::ofstream logfile("grid_tools_02/output");



template <int dim>
void test1 ()
{
				   // test 1: hypercube
  if (true)
    {
      Triangulation<dim,dim+1> tria;
      GridGenerator::hyper_cube(tria);

      for (unsigned int i=0; i<2; ++i)
	{
	  tria.refine_global(2);
	  deallog << dim << "d, "
		  << "hypercube diameter, "
		  << i*2
		  << " refinements: "
		  << GridTools::diameter (tria)
		  << std::endl;
	}
    }
}



int main ()
{
  deallog << std::setprecision(4);
  logfile << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test1<1> ();
  test1<2> ();

  return 0;
}

