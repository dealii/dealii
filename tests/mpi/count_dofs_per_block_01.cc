// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// test DoFTools::count_dofs_per_block in parallel


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>
#include <cstdlib>
#include <numeric>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim>
  triangulation (MPI_COMM_WORLD,
                 Triangulation<dim>::limit_level_difference_at_vertices);

  FESystem<dim> fe (FE_Q<dim>(3),2,
                    FE_DGQ<dim>(1),1);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (4);

  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  std::vector<types::global_dof_index> dofs_per_block (fe.n_blocks());
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    for (unsigned int i=0; i<fe.n_blocks(); ++i)
      deallog << "Block " << i << " has " << dofs_per_block[i] << " dofs" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }

}
