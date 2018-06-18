// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// check LA::Vector::compress(VectorOperation::min/max) from ghosts

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;


  // each processor owns 2 indices and all
  // are ghosting element 1 (the second)
  IndexSet local_owned(numproc * 2);
  local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);

  // create vector
  LinearAlgebra::distributed::Vector<double> v(local_owned,
                                               local_relevant,
                                               MPI_COMM_WORLD);

  // set local values
  v(myid * 2)     = myid * 2.0;
  v(myid * 2 + 1) = myid * 2.0 + 1.0;
  v.compress(VectorOperation::add);
  v *= 2.0;

  // check setup of vectors
  deallog << myid << ":"
          << "first owned entry: " << v(myid * 2) << std::endl;
  deallog << myid << ":"
          << "second owned entry: " << v(myid * 2 + 1) << std::endl;

  // set ghost dof on not owning processor and maximize
  if (myid)
    v(1) = 7. * myid;
  v.compress(VectorOperation::max);

  // import ghosts onto all procs
  v.update_ghost_values();

  // check
  deallog << myid << ":"
          << "ghost entry: " << v(1) << std::endl;

  // ghosts are set to zero
  v.zero_out_ghosts();

  // minimize
  v.compress(VectorOperation::min);
  v.update_ghost_values();

  // check
  deallog << myid << ":"
          << "ghost entry: " << v(1) << std::endl;

  // update of ghost value from owner and minimize
  v.zero_out_ghosts();
  if (!myid)
    v(1) = -1.;
  v.compress(VectorOperation::min);
  v.update_ghost_values();

  // check
  deallog << myid << ":"
          << "ghost entry: " << v(1) << std::endl;

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;

  test();
}
