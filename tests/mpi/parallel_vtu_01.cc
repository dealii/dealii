// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// create a parallel DoFHandler and output data using the parallel vtk output
// (uses MPI IO)

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(2);
  dofh.distribute_dofs (fe);


  TrilinosWrappers::MPI::Vector x;
  x.reinit(dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  x=2.0;

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dofh);
  data_out.add_data_vector (x, "x");
  data_out.build_patches ();

  data_out.write_vtu_in_parallel("output.vtu", MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (myid==0)
    {
      std::vector<std::string> filenames;
      filenames.push_back("output.vtu");
      {	
	std::ofstream master("output.pvtu");
	data_out.write_pvtu_record (master, filenames);
      }
      
      cat_file("output.vtu");
      cat_file("output.pvtu");      
    }
  deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  test<2>();
}
