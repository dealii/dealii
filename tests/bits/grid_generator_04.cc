//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test GridGenerator::hyper_cube_with_hole

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


template<int dim>
void test(std::ostream& out)
{
  GridOut go;
  go.set_flags(GridOutFlags::Ucd(false, true));
  Triangulation<dim> tr;
  
  std::vector<double> radii;
  radii.push_back(.2);
  radii.push_back(.3);
  radii.push_back(.4);
  
  
  std::vector<double> radiiext;
  radiiext.push_back(.3);
  radiiext.push_back(.4);
  radiiext.push_back(.5);
  
  std::vector<double> Ls;
  Ls.push_back(.1);
  Ls.push_back(1);
  Ls.push_back(10);

  std::vector<unsigned int> Nzs;
  Nzs.push_back(2);
  Nzs.push_back(3);
  Nzs.push_back(4);

  for(unsigned int i=0; i<radii.size(); ++i)
    for(unsigned int k=0; k< (dim == 2 ? 1 : Ls.size()); ++k)
      for(unsigned int l=0; l< (dim == 2 ? 1 : Ls.size()); ++l) {
	  
	  out << "               ====================" << std::endl
	      << "Inner radius = " << radii[i] << std::endl
	      << "Outer radius = " << radiiext[i] << std::endl
	      << "Depth        = " << Ls[k] << std::endl
	      << "Nzs	       = " << Nzs[l] << std::endl
	      << "No colorize    ====================" << std::endl;

	  // No colorize
	  try {
	    GridGenerator::hyper_cube_with_cylindrical_hole(tr, radii[i], radiiext[i], Ls[k], Nzs[l], false);
	  } catch(...) {
	    out << "Generation failed." << std::endl;
	  }
	  if (tr.n_cells() > 0)
	    go.write_ucd(tr, out);
	  
	  tr.clear();
	  
	  out << "Colorize       ====================" << std::endl;
	  try {
	    GridGenerator::hyper_cube_with_cylindrical_hole(tr, radii[i], radiiext[i], Ls[k], Nzs[l], true);
	  } catch(...) {
	    out << "Generation failed." << std::endl;
	  }
	  if (tr.n_cells() > 0)
	    go.write_ucd(tr, out);
	  tr.clear();
	  
	  out << "               ====================" << std::endl;
	}
}


int main()
{
  std::ofstream logfile("grid_generator_04/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
