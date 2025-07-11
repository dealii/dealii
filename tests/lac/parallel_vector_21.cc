// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that the data range representing ghosts is really initialized to zero
// when doing reinit() from another vector and manually setting the local
// range

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
  unsigned int my_id   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  IndexSet locally_owned(n_procs * 2);
  locally_owned.add_range(my_id * 2, my_id * 2 + 2);
  IndexSet ghost_set(n_procs * 2);
  ghost_set.add_index(0);
  ghost_set.add_index(2);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    locally_owned, ghost_set, MPI_COMM_WORLD);

  // create vector without actually setting the entries since they will be
  // overwritten soon anyway
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v2;
  v2.reinit(v, true);

  // set locally owned range of v2 manually
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> v2_view(
    v2.get_values(), v2.locally_owned_size());
  Kokkos::deep_copy(v2_view, 1.);

  // add entries to ghost values
  // Because of limitation in import, the IndexSet of the ReadWriteVector needs
  // to have the local elements.
  IndexSet workaround_set(locally_owned);
  workaround_set.add_index(0);
  workaround_set.add_index(2);
  workaround_set.compress();
  LinearAlgebra::ReadWriteVector<double> rw_vector(workaround_set);
  rw_vector(0) += 1.;
  rw_vector(2) += 1.;
  v2.import_elements(rw_vector, VectorOperation::add);

  // now we should have the correct data, not some uninitialized trash that
  // resided in the ghost range
  v2.print(deallog.get_file_stream());

  v2.update_ghost_values();
  v2.print(deallog.get_file_stream());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;
  test();
}
