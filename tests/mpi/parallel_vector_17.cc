// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
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


// check LinearAlgebra::distributed::Vector::partitioners_are_compatible and
// partitioners_are_globally_compatible

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // all processors own 2 elements
  IndexSet local_owned(numproc*2);
  local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(1,2);

  LinearAlgebra::distributed::Vector<double> v1, v2, v3, v4, v5, v6;
  v1.reinit(local_owned, MPI_COMM_WORLD);
  v2.reinit(local_owned, local_relevant, MPI_COMM_WORLD);
  v3.reinit(local_owned, local_relevant, MPI_COMM_WORLD);
  v4.reinit(v3);
  IndexSet local_owned_5(numproc*3);
  local_owned_5.add_range(myid*3,myid*3+3);
  v5.reinit(local_owned_5, MPI_COMM_WORLD);

  // create sub-communicator on two processors
  MPI_Comm newcomm;
  MPI_Comm_split(MPI_COMM_WORLD, myid/2, myid%2, &newcomm);
  IndexSet local_owned_b(6);
  if (myid<2)
    local_owned_b.add_range(myid*3,myid*3+3);
  v6.reinit(local_owned_b, newcomm);

  deallog << "local compatibility  v1-v2: "
          << v1.partitioners_are_compatible(*v2.get_partitioner()) << " "
          << v2.partitioners_are_compatible(*v1.get_partitioner()) << std::endl;
  deallog << "global compatibility v1-v2: "
          << v1.partitioners_are_globally_compatible(*v2.get_partitioner()) << " "
          << v2.partitioners_are_globally_compatible(*v1.get_partitioner()) << std::endl;
  deallog << "local compatibility  v2-v3: "
          << v2.partitioners_are_compatible(*v3.get_partitioner()) << " "
          << v3.partitioners_are_compatible(*v2.get_partitioner()) << std::endl;
  deallog << "global compatibility v2-v3: "
          << v2.partitioners_are_globally_compatible(*v3.get_partitioner()) << " "
          << v3.partitioners_are_globally_compatible(*v2.get_partitioner()) << std::endl;
  deallog << "local compatibility  v3-v4: "
          << v4.partitioners_are_compatible(*v3.get_partitioner()) << " "
          << v3.partitioners_are_compatible(*v4.get_partitioner()) << std::endl;
  deallog << "global compatibility v3-v4: "
          << v4.partitioners_are_globally_compatible(*v3.get_partitioner()) << " "
          << v3.partitioners_are_globally_compatible(*v4.get_partitioner()) << std::endl;
  deallog << "local compatibility  v4-v5: "
          << v4.partitioners_are_compatible(*v5.get_partitioner()) << " "
          << v5.partitioners_are_compatible(*v4.get_partitioner()) << std::endl;
  deallog << "global compatibility v4-v5: "
          << v4.partitioners_are_globally_compatible(*v5.get_partitioner()) << " "
          << v5.partitioners_are_globally_compatible(*v4.get_partitioner()) << std::endl;
  deallog << "local compatibility  v5-v6: "
          << v6.partitioners_are_compatible(*v5.get_partitioner()) << " "
          << v5.partitioners_are_compatible(*v6.get_partitioner()) << std::endl;
  deallog << "global compatibility v5-v6: "
          << v6.partitioners_are_globally_compatible(*v5.get_partitioner()) << " "
          << v5.partitioners_are_globally_compatible(*v6.get_partitioner()) << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());
  MPILogInitAll log;

  test();

}
