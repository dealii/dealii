// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check that indeed Triangulation::create_triangulation throws an
// exception if we have distorted cells. this test is like the _01
// case except that it reads the mesh from a file

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>


template <int dim>
void check ()
{
  Triangulation<dim> coarse_grid (Triangulation<dim>::none, true);

  GridIn<dim> gi;
  std::ifstream in (dim == 2 ? SOURCE_DIR "/grids/2d" :
      SOURCE_DIR "/grids/3d");

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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}



