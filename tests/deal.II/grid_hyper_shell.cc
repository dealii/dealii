//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


#include "../tests.h"
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>
#include <base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <iomanip>

std::ofstream logfile("grid_hyper_shell/output");


template<int dim>
void check (double r1, double r2, unsigned int n, bool)
{
  Point<dim> center;
  Triangulation<dim> tria;
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
	  Assert (subset.distorted_cells.size() == 0,
		  ExcInternalError());
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
  check<2> (4., 5., 10, true);
  check<3> (3., 5., 6, true);
}
