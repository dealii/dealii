// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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

std::ofstream logfile("output");


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
