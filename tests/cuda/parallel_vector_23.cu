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
#include <deal.II/lac/read_write_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


__global__ void
set_value(double *values_dev, unsigned int index, double val)
{
  values_dev[index] = val;
}


template <typename Number>
double
print_value(Number *values_dev, unsigned int index)
{
  static std::vector<Number> cpu_value(1);
  Utilities::CUDA::copy_to_host(values_dev + index, cpu_value);
  return cpu_value[0];
}



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
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> v(
    local_owned, local_relevant, MPI_COMM_WORLD);
  const auto &partitioner = v.get_partitioner();

  // the read write vector additionally has ghost elements
  IndexSet                               read_write_owned(numproc * 2);
  LinearAlgebra::ReadWriteVector<double> read_write_vector(local_relevant);

  read_write_vector.local_element(0) = myid;
  read_write_vector(1)               = 2. * myid;

  v.import(read_write_vector, VectorOperation::max);
  v.update_ghost_values();

  deallog << myid << ":"
          << "ghost entry after max: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  if (myid == 0)
    read_write_vector(1) = -1.0;

  v.import(read_write_vector, VectorOperation::min);
  v.update_ghost_values();

  deallog << myid << ":"
          << "ghost entry after min: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;


  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;

  init_cuda(true);

  test();
}
