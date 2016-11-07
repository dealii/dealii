// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// check memory consumption parallel vector

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // each processor from processor 1 to 8 owns 2 indices (the other processors
  // do not own any dof), and all processors are ghosting element 1
  IndexSet local_owned(30);
  if (myid == 0)
    local_owned.add_range(0,20);
  else
    local_owned.add_range(std::min(20U,20*myid),30U);
  IndexSet local_relevant(30);
  local_relevant = local_owned;
  local_relevant.add_range(19,20);

  LinearAlgebra::distributed::Vector<double> v(local_owned, local_relevant,
                                               MPI_COMM_WORLD);

  deallog << v.memory_consumption() << std::endl;

  LinearAlgebra::distributed::BlockVector<double> w(3);
  for (unsigned int i=0; i<3; ++i)
    {
      w.block(i) = v;
      w.block(i) *= (i+1);
    }
  w.collect_sizes();

  deallog << w.memory_consumption() << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  MPILogInitAll log;
  test();
}
