// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <fstream>


template<int dim>
void test()
{
  // create 3 triangulations
  Triangulation<dim> tria[3];

  GridGenerator::hyper_cube (tria[0]);
  tria[0].refine_global (1);

  GridGenerator::hyper_cube (tria[1]);
  GridTools::scale (2, tria[1]);
  tria[1].refine_global (2);

  if (dim != 1)
    GridGenerator::hyper_ball (tria[2]);
  else
    {
      GridGenerator::hyper_cube (tria[2]);
      GridTools::shift (Point<dim>(2.), tria[2]);
    }

  tria[2].refine_global (3);

  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=0; j<3; ++j)
      {
        Assert (GridTools::have_same_coarse_mesh (tria[i], tria[j])
                ==
                (i == j),
                ExcInternalError());

        deallog << "meshes " << i << " and " << j << ": "
                << (GridTools::have_same_coarse_mesh (tria[i], tria[j])
                    ?
                    "true"
                    :
                    "false")
                << std::endl;
      }
}


int main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

