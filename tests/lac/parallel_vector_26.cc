// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check for a bug where compress cannot be called for a LA::d::Vector<Host> in
// a Device file

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.templates.h>

#include "../tests.h"


void
test()
{
  unsigned int       rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int total_size = 100;
  const unsigned int ghost_size = 75;
  const unsigned int local_size = 50;

  IndexSet local_owned(total_size);
  if (rank == 0)
    local_owned.add_range(0, local_size);
  else
    local_owned.add_range(total_size - local_size, total_size);

  IndexSet ghost_indices(total_size);
  if (rank == 0)
    ghost_indices.add_range(0, ghost_size);
  else
    ghost_indices.add_range(total_size - ghost_size, total_size);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Host> v(
    local_owned, ghost_indices, MPI_COMM_WORLD);

  for (unsigned int i = 0; i < ghost_size; ++i)
    v.local_element(i) = i;
  v.compress(VectorOperation::add);


  if (rank == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(rank));

  if (rank == 0)
    {
      initlog();
      deallog << std::setprecision(4);
      test();
    }
  else
    test();
}
