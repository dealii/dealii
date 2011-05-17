//----------------------------  grid_in_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in_02.cc  ---------------------------

// when reading a particular input file with a mesh, it turns out that
// there are two cells where we list the same cell as neighbor through two
// different faces. this of course can't be
//
// the actual cause turned out to be that in each of these 2 cases, the
// (machine generated) input file had two cells that shared 3 (!) vertices.
// in each case, one of the two cells has been removed from the input
// file to fix the testcase.

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

std::ofstream logfile("grid_in_02/output");


template <int dim>
void test2 ()
{
                                   // read a much larger grid (30k
                                   // cells). with the old grid
                                   // reordering scheme, this took >90
                                   // minutes (exact timing not
                                   // available, program was killed
                                   // before), with the new one it
                                   // takes less than 8 seconds
  Triangulation<dim> tria (Triangulation<dim>::none, true);
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in ("grid_in_02/2d.xda");
  try
    {
      gi.read_xda (in);
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
				       // ignore the exception that we
				       // get because the mesh has
				       // distorted cells
      deallog << dcv.distorted_cells.size() << " cells are distorted."
	      << std::endl;
    }

  Triangulation<2>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();
  for (; cell != endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      for (unsigned int e=0; e<GeometryInfo<2>::faces_per_cell; ++e)
	if (f != e)
	  if (!cell->at_boundary(e) && !cell->at_boundary(f))
	    Assert (cell->neighbor(e) !=
		    cell->neighbor(f),
		    ExcInternalError());
  deallog << "OK" << std::endl;
}




int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test2<2> ();
}

