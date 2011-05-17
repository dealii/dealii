//----------------------------  have_same_coarse_mesh_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  have_same_coarse_mesh_01.cc  ---------------------------
// check GridTools::have_same_coarse_mesh for triangulations


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
  std::ofstream logfile ("have_same_coarse_mesh_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

