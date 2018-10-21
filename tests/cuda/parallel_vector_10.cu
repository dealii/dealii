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


// check operator= when we do some operations with ghosts

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


  // each processor owns 2 indices and all
  // are ghosting element 1 (the second)
  IndexSet local_owned(numproc * 2);
  local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);

  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> v(
    local_owned, local_relevant, MPI_COMM_WORLD);
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> w(v);

  // set local values and check them
  LinearAlgebra::ReadWriteVector<double> rw_vector(local_owned);
  rw_vector(myid * 2)     = myid * 2.0;
  rw_vector(myid * 2 + 1) = myid * 2.0 + 1.0;
  v.import(rw_vector, VectorOperation::insert);

  v.update_ghost_values();

  // check that the value of the ghost is 1.0
  IndexSet ghost_set(numproc * 2);
  ghost_set.add_index(1);
  LinearAlgebra::ReadWriteVector<double> ghost_vector(ghost_set);
  ghost_vector.import(v, VectorOperation::insert);
  DEAL_II_AssertThrow(ghost_vector(1) == 1., ExcInternalError());

  // copy vector
  w = v;
  v *= 2.0;

  v.update_ghost_values();
  w.update_ghost_values();
  ghost_vector.import(v, VectorOperation::insert);
  DEAL_II_AssertThrow(ghost_vector(1) == 2., ExcInternalError());
  ghost_vector.import(w, VectorOperation::insert);
  DEAL_II_AssertThrow(ghost_vector(1) == 1., ExcInternalError());

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

  Utilities::CUDA::Handle cuda_handle;
  // By default, all the ranks will try to access the device 0. This is fine if
  // we have one rank per node _and_ one gpu per node. If we have multiple GPUs
  // on one node, we need each process to access a different GPU. We assume that
  // each node has the same number of GPUs.
  int         n_devices       = 0;
  cudaError_t cuda_error_code = cudaGetDeviceCount(&n_devices);
  DEAL_II_AssertCuda(cuda_error_code);
  int device_id   = myid % n_devices;
  cuda_error_code = cudaSetDevice(device_id);
  DEAL_II_AssertCuda(cuda_error_code);

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
