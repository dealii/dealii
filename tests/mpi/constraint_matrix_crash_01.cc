// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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



// ConstraintMatrix::add_line crashes in release mode (missing compress inside ConstraintMatrix)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/parallel_block_vector.h>

#include <fstream>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  IndexSet local_active_together(3);
  local_active_together.add_range(0,3);
  //local_active_together.compress();
  
  ConstraintMatrix cm(local_active_together);
  cm.add_line(1);
  cm.close();
  deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  test();
  return 0;
}
