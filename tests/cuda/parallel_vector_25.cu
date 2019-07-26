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


// check import from LA::d::Vector to another LA::d::Vector

#include <deal.II/base/cuda.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

void
test()
{
  unsigned int       rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int size = 100;
  const unsigned int local_size = 50;
  const unsigned int ghost_size = 75;

  IndexSet locally_owned(size);
  if (rank == 0)
    locally_owned.add_range(0, local_size);
  else
    locally_owned.add_range(size - local_size, size);
  locally_owned.compress();

  IndexSet ghost_indices(size);
  if (rank == 0)
    ghost_indices.add_range(0, ghost_size);
  else
    ghost_indices.add_range(size - ghost_size, size);
  ghost_indices.size();

  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> vec_ref(
    locally_owned, ghost_indices, MPI_COMM_WORLD);
  for (unsigned int i = 0; i < local_size; ++i)
    vec_ref.local_element(i) = i;
  vec_ref.compress(VectorOperation::insert);

  auto partitioner = vec_ref.get_partitioner();
  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> vec_dev(
    partitioner);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> vec_host(
    partitioner);

  // Assignment from Host to CUDA
  vec_dev.import(vec_ref, VectorOperation::insert);

  // Assignment from CUDA to Host
  vec_host.import(vec_dev, VectorOperation::insert);

  for (unsigned int i = 0; i < ghost_size; ++i)
    {
      AssertThrow(std::fabs(vec_ref.local_element(i) -
                            vec_host.local_element(i)) < 1e-12,
                  ExcInternalError());
    }

  if (rank == 0)
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

      test();
    }
  else
    test();
}
