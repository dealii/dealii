// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check LinearAlgebra::distributed::Vector's move constructor. Should work like
// swap.

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

DeclException2(ExcNonEqual,
               double,
               double,
               << "Left compare: " << arg1 << ", right compare: " << arg2);

void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // vector 0:
  // global size: 20, locally_owned_size: 3 as long as
  // less than 20
  const unsigned int locally_owned_size0 = 3;
  const unsigned int global_size0 =
    std::min(20U, locally_owned_size0 * numproc);
  const unsigned int my_start0 =
    std::min(locally_owned_size0 * myid, global_size0);
  const unsigned int my_end0 =
    std::min(locally_owned_size0 * (myid + 1), global_size0);
  const unsigned int actual_locally_owned_size0 = my_end0 - my_start0;

  IndexSet local_owned0(global_size0);
  if (my_end0 > my_start0)
    local_owned0.add_range(static_cast<unsigned int>(my_start0),
                           static_cast<unsigned int>(my_end0));
  IndexSet local_relevant0(global_size0);
  local_relevant0 = local_owned0;
  local_relevant0.add_index(2);
  if (numproc > 2)
    local_relevant0.add_index(8);

  LinearAlgebra::distributed::Vector<double> v0(local_owned0,
                                                local_relevant0,
                                                MPI_COMM_WORLD);
  v0 = 1;
  // check assignment in initial state
  for (unsigned int i = 0; i < v0.locally_owned_size(); ++i)
    AssertThrow(v0.local_element(i) == 1.,
                ExcNonEqual(v0.local_element(i), 1.));

  // check ghost elements in initial state
  v0.update_ghost_values();
  AssertThrow(v0(2) == 1., ExcNonEqual(v0(2), 1.));
  if (numproc > 2)
    AssertThrow(v0(8) == 1., ExcNonEqual(v0(8), 1.));
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Initial set and ghost update OK" << std::endl;

  // now move.
  LinearAlgebra::distributed::Vector<double> v1 = std::move(v0);
  // Make sure we actually moved and not copied
  AssertDimension(v0.locally_owned_size(), 0);
  AssertDimension(v1.locally_owned_size(), actual_locally_owned_size0);
  AssertDimension(v0.size(), 0);
  AssertDimension(v1.size(), global_size0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "First move: dimensions OK" << std::endl;
  for (unsigned int i = 0; i < actual_locally_owned_size0; ++i)
    AssertThrow(v1.local_element(i) == 1.,
                ExcNonEqual(v1.local_element(i), 1.));
  // Since we moved the ghost values should be present
  for (const auto &ghost_index : v1.get_partitioner()->ghost_indices())
    {
      AssertThrow(v1(ghost_index) == 1., ExcNonEqual(v1(ghost_index), 1.));
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "First move: local values OK" << std::endl;
  v0.update_ghost_values();
  v1.update_ghost_values();
  AssertThrow(v1(2) == 1., ExcNonEqual(v1(2), 1.));
  if (numproc > 2)
    AssertThrow(v1(8) == 1., ExcNonEqual(v1(8), 1.));
  if (numproc > 2)
    AssertThrow(v1(8) == 1., ExcNonEqual(v1(8), 1.));
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Ghost values after first move OK" << std::endl;
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
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
