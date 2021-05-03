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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// check that handling of ghost elements in parallel distributed vectors works
// appropriately when creating a vector from a non-ghosted source vector using
// the assignment operator

#include <deal.II/base/cuda.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

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

  // processor 0 and 1 own 2 indices each, higher processors nothing, all are
  // ghosting global elements 1 and 3
  IndexSet local_owned(std::min(numproc * 2, 4U));
  if (myid < 2)
    local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);
  if (numproc > 1)
    local_relevant.add_range(3, 4);

  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> v(
    local_owned, local_relevant, MPI_COMM_WORLD);

  // set local values
  LinearAlgebra::ReadWriteVector<double> rw_vector(local_owned);
  if (myid < 2)
    {
      rw_vector(myid * 2)     = myid * 2.0;
      rw_vector(myid * 2 + 1) = myid * 2.0 + 1.0;
    }

  v.import(rw_vector, VectorOperation::insert);

  if (myid == 0)
    deallog << "v has ghost elements: " << v.has_ghost_elements() << std::endl;

  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> w, x;
  w = v;
  if (myid == 0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  v.update_ghost_values();
  w = v;
  if (myid == 0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  v.zero_out_ghost_values();
  w = v;
  if (myid == 0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  w.zero_out_ghost_values();
  w = v;
  if (myid == 0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  v.update_ghost_values();
  x = v;
  if (myid == 0)
    deallog << "x has ghost elements: " << x.has_ghost_elements() << std::endl;

  x.zero_out_ghost_values();
  if (myid == 0)
    deallog << "x has ghost elements: " << x.has_ghost_elements() << std::endl;

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  init_cuda(true);

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
