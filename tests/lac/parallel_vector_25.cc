// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check import from LA::d::Vector to another LA::d::Vector

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
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> vec_dev(
    partitioner);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> vec_host(
    partitioner);

  // Assignment from Host to Default
  vec_dev.import_elements(vec_ref, VectorOperation::insert);

  // Assignment from Default to Host
  vec_host.import_elements(vec_dev, VectorOperation::insert);

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

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
