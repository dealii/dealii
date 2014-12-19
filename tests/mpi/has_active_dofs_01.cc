// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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



// verify that DoFHandler::has_active_dofs() works also on meshes
// where some processors have no active dofs because they don't own
// any cells

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>


#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  // create a mesh with fewer cells than there are MPI processes
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  // assign DoFs to the mesh
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tr);
  dof_handler.distribute_dofs (fe);
  
  // there are 2^(2*dim) cells to be owned by this distributed
  // triangulation. however, since p4est makes sure that all active
  // children belonging to the same parent are stored on the same
  // processor, only 2^dim processors can own cells.
  //
  // let each of the processors report whether they own
  // cells or now by setting a bit mask that we then add up.
  // output which processor owns something and which don't
  const int cells_owned = Utilities::MPI::sum (tr.n_locally_owned_active_cells() > 0
					       ?
					       1 << Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)
					       :
					       0,
					       MPI_COMM_WORLD);
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
	deallog << "Processor " << i << " has cells: "
		<< ((cells_owned & (1 << i)) ? "yes" : "no")
		<< std::endl;
      
      unsigned int n_owning_processors = 0;      
      for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
	if (cells_owned & (1 << i))
	  ++n_owning_processors;
      deallog << "# processors owning cells=" << n_owning_processors << std::endl;
    }

  // now check how many processors own DoFs using the same procedure
  const int dofs_owned = Utilities::MPI::sum (dof_handler.has_active_dofs()
					      ?
					      1 << Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)
					      :
					      0,
					      MPI_COMM_WORLD);
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
	deallog << "Processor " << i << " has dofs: "
		<< ((dofs_owned & (1 << i)) ? "yes" : "no")
		<< std::endl;
      
      unsigned int n_owning_processors = 0;      
      for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
	if (dofs_owned & (1 << i))
	  ++n_owning_processors;
      deallog << "# processors owning dofs=" << n_owning_processors << std::endl;
    }
  
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
