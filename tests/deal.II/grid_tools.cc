//----------------------------  grid_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_test.cc  ---------------------------


#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>

#include <fstream>


std::ofstream logfile("grid_tools.output");



template <int dim>
void test ()
{
				   // test 1: hypercube
  if (true)
    {
      Triangulation<dim> tria;
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
	};
    };

				   // test 2: hyperball
  if (dim == 2)
    {
      Triangulation<dim> tria;
      GridGenerator::hyper_ball(tria, Point<dim>(), 1);

      for (unsigned int i=0; i<2; ++i)
	{
	  tria.refine_global(2);
	  deallog << dim << "d, "
		  << "hyperball diameter, "
		  << i*2
		  << " refinements: "
		  << GridTools::diameter (tria)
		  << std::endl;
	};
    };
};


int main ()
{
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
  
  return 0;
};
