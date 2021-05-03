// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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


// build a vector whose elements exceed the size of unsigned int in case of 64
// bit indices. To avoid excessive memory consumption, let the vector start at
// a number close to the maximum of unsigned int but extend past the last
// index

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


__global__ void
set_value(double *values_dev, unsigned int index, double val)
{
  values_dev[index] = val;
}

void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  types::global_dof_index min_index  = 0xffffffffU - 39;
  types::global_dof_index local_size = 42;
  IndexSet                local_owned(min_index + numproc * local_size);
  local_owned.add_range(min_index + myid * local_size,
                        min_index + (myid + 1) * local_size);

  // all processors ghost some entries around invalid_unsigned_int and on the
  // border between two processors
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(min_index + 38, min_index + 40);
  local_relevant.add_range(min_index + 41, min_index + 43);

  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> v(
    local_owned, local_relevant, MPI_COMM_WORLD);

  deallog << "Local range of proc 0: " << v.local_range().first << " "
          << v.local_range().second << std::endl;

  // set local values
  for (types::global_dof_index i = 0; i < local_size; ++i)
    {
      double *values_dev = v.get_values();
      set_value<<<1, 1>>>(values_dev, i, min_index + myid * local_size + i);
    }

  deallog << "vector norm: " << v.l2_norm() << std::endl;

  // check ghost values
  v.print(deallog.get_file_stream(), 12, false, false);
  v.update_ghost_values();
  v.print(deallog.get_file_stream(), 12, false, false);

  v.zero_out_ghost_values();
  double *    values_dev  = v.get_values();
  const auto &partitioner = v.get_partitioner();
  set_value<<<1, 1>>>(values_dev,
                      partitioner->global_to_local(min_index + 38),
                      min_index);
  set_value<<<1, 1>>>(values_dev,
                      partitioner->global_to_local(min_index + 39),
                      min_index * 2);
  set_value<<<1, 1>>>(values_dev,
                      partitioner->global_to_local(min_index + 41),
                      min_index + 7);
  set_value<<<1, 1>>>(values_dev,
                      partitioner->global_to_local(min_index + 42),
                      -static_cast<double>(min_index));
  v.compress(VectorOperation::add);
  v.update_ghost_values();
  v.print(deallog.get_file_stream(), 12, false, false);

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll mpilog;

  init_cuda(true);

  deallog << std::setprecision(12);

  test();
}
