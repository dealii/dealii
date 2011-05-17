//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// verify the distortion in cells of a hyper shell with 6 cells upon
// refinement

#include "../tests.h"
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

std::ofstream logfile("grid_hyper_shell/output");


template<int dim>
void check (double r1, double r2, unsigned int n)
{
  deallog << "dim=" << dim << std::endl;
  
  Point<dim> center;
  Triangulation<dim> tria (Triangulation<dim>::none, true);
  GridGenerator::hyper_shell (tria, center, r1, r2, n, true);
  static const HyperShellBoundary<dim> boundary(center);
  tria.set_boundary(0, boundary);

  for (unsigned int i=0; i<2; ++i)
    {
      try
	{
	  tria.refine_global(1);
	}
      catch (typename Triangulation<dim>::DistortedCellList &dcv)
	{
	  deallog << "Found " << dcv.distorted_cells.size()
		  << " distorted cells" << std::endl;

	  typename Triangulation<dim>::DistortedCellList
	    subset = GridTools::fix_up_distorted_child_cells (dcv,
							      tria);
	  deallog << subset.distorted_cells.size()
		  << " distorted cells remaining" << std::endl;
	}
    }

  GridOut grid_out;
  GridOutFlags::DX flags;
  flags.write_faces = true;
  grid_out.set_flags(flags);
  grid_out.write_dx (tria, logfile);
}


int main()
{
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  check<2> (4., 5., 10);
  check<3> (3., 5., 6);
}
