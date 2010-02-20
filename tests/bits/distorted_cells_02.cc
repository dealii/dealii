//----------------------------  distorted_cells_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  distorted_cells_02.cc  ---------------------------


// check that indeed Triangulation::create_triangulation throws an
// exception if we have distorted cells. this test is like the _01
// case except that it reads the mesh from a file

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_in.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <fstream>


template <int dim>
void check ()
{
  Triangulation<dim> coarse_grid (Triangulation<dim>::none, true);

  GridIn<dim> gi;
  std::ifstream in (dim == 2 ? "distorted_cells_02/2d" : "distorted_cells_02/3d");

  gi.attach_triangulation (coarse_grid);


  bool flag = false;
  try
    {
      gi.read_ucd (in);
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      flag = true;

      deallog << dcv.distorted_cells.size() << " distorted cells"
	      << std::endl;
      Assert (dcv.distorted_cells.front() == coarse_grid.begin(0),
	      ExcInternalError());
    }

  Assert (flag == true, ExcInternalError());
  Assert (coarse_grid.n_levels() == 1, ExcInternalError());
  Assert (coarse_grid.n_active_cells() == 1, ExcInternalError());
}


int main ()
{
  std::ofstream logfile("distorted_cells_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}



