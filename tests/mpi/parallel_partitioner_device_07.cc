// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for the Partitioner with a smaller ghost index set within a larger one
// regarding the import_from_ghosted_array() calls
#include <deal.II/base/index_set.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/utilities.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../tests.h"

template <typename Number>
void
print_device_view(
  const Kokkos::View<Number *, MemorySpace::Default::kokkos_space> device_view)
{
  std::vector<Number> cpu_values(device_view.size());
  Kokkos::deep_copy(Kokkos::View<Number *, Kokkos::HostSpace>(
                      cpu_values.data(), cpu_values.size()),
                    device_view);
  for (Number value : cpu_values)
    deallog << value << " ";
  deallog << std::endl;
}


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  Assert(numproc > 2, ExcNotImplemented());

  const unsigned int set = 50;
  AssertIndexRange(numproc, set - 2);
  const unsigned int      local_size  = set - myid;
  types::global_dof_index global_size = 0;
  types::global_dof_index my_start    = 0;
  for (unsigned int i = 0; i < numproc; ++i)
    {
      global_size += set - i;
      if (i < myid)
        my_start += set - i;
    }

  // each processor owns some indices and all are ghosting elements from three
  // processors (the second). some entries are right around the border between
  // two processors
  IndexSet local_owned(global_size);
  local_owned.add_range(my_start, my_start + local_size);
  IndexSet local_relevant_1(global_size), local_relevant_2(global_size);
  local_relevant_1                          = local_owned;
  types::global_dof_index ghost_indices[10] = {1,
                                               2,
                                               13,
                                               set - 2,
                                               set - 1,
                                               set,
                                               set + 1,
                                               2 * set,
                                               2 * set + 1,
                                               2 * set + 3};
  local_relevant_1.add_indices(&ghost_indices[0], ghost_indices + 10);
  if (myid > 0)
    local_relevant_1.add_range(my_start - 10, my_start);
  if (myid < numproc - 1)
    local_relevant_1.add_range(my_start + local_size,
                               my_start + local_size + 10);

  local_relevant_2 = local_owned;
  local_relevant_2.add_indices(&ghost_indices[0], ghost_indices + 10);
  if (myid > 0)
    local_relevant_2.add_index(my_start - 10);
  if (myid < numproc - 1)
    local_relevant_2.add_index(my_start + local_size + 9);

  Utilities::MPI::Partitioner v(local_owned, local_relevant_1, MPI_COMM_WORLD);
  Utilities::MPI::Partitioner w(local_owned, MPI_COMM_WORLD);
  w.set_ghost_indices(local_relevant_2, v.ghost_indices());

  IndexSet local_relevant_3(global_size);
  local_relevant_3.add_index(2);
  if (myid > 0 && my_start > 0)
    local_relevant_3.add_range(my_start - 10, my_start);
  Utilities::MPI::Partitioner x(local_owned, MPI_COMM_WORLD);
  x.set_ghost_indices(local_relevant_3, v.ghost_indices());

  // set up a ghost array with some entries
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghost_array(
    "ghost_array", v.n_ghost_indices());
  ArrayView<unsigned int, MemorySpace::Default> ghost_array_view(
    ghost_array.data(), ghost_array.size());
  Kokkos::deep_copy(ghost_array, 1);

  // set up other arrays
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space>
    locally_owned_array("locally_owned_array", local_size);
  ArrayView<unsigned int, MemorySpace::Default> locally_owned_array_view(
    locally_owned_array.data(), locally_owned_array.size());

  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> temp_array(
    "temp_array", v.n_import_indices());
  ArrayView<unsigned int, MemorySpace::Default> temp_array_view(
    temp_array.data(), temp_array.size());

  std::vector<MPI_Request> requests;

  // send the full array
  {
    Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghosts(
      "ghosts", ghost_array_view.size());
    ArrayView<unsigned int, MemorySpace::Default> ghosts_view(ghosts.data(),
                                                              ghosts.size());
    Kokkos::deep_copy(ghosts, ghost_array);

    v.import_from_ghosted_array_start<unsigned int, MemorySpace::Default>(
      VectorOperation::add, 3, ghosts_view, temp_array_view, requests);
    v.import_from_ghosted_array_finish<unsigned int, MemorySpace::Default>(
      VectorOperation::add,
      temp_array_view,
      locally_owned_array_view,
      ghosts_view,
      requests);
    // check that the ghost entries are zeroed out in these calls
    deallog << "v ghost entries (should be zero up to index "
            << v.n_ghost_indices() - 1 << "):" << std::endl;
    print_device_view(ghosts);
  }
  deallog << "From all ghosts: ";
  print_device_view(locally_owned_array);

  // send only the array in w
  Kokkos::deep_copy(locally_owned_array, 0);
  Assert(temp_array_view.size() >= w.n_import_indices(), ExcInternalError());
  ArrayView<unsigned int, MemorySpace::Default> temp_array_view_w(
    temp_array_view.data(), w.n_import_indices());
  {
    Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghosts(
      "ghosts", ghost_array_view.size());
    ArrayView<unsigned int, MemorySpace::Default> ghosts_view(ghosts.data(),
                                                              ghosts.size());
    Kokkos::deep_copy(ghosts, ghost_array);

    w.import_from_ghosted_array_start<unsigned int, MemorySpace::Default>(
      VectorOperation::add, 3, ghosts_view, temp_array_view_w, requests);
    w.import_from_ghosted_array_finish<unsigned int, MemorySpace::Default>(
      VectorOperation::add,
      temp_array_view_w,
      locally_owned_array_view,
      ghosts_view,
      requests);

    // check that the ghost entries are zeroed out in these calls
    deallog << "w ghost entries (should be zero up to index "
            << w.n_ghost_indices() - 1 << "):" << std::endl;
    print_device_view(ghosts);
  }
  deallog << "From reduced ghosts 1: ";
  print_device_view(locally_owned_array);

  // send only the array in x
  Kokkos::deep_copy(locally_owned_array, 0);
  Assert(temp_array_view.size() >= x.n_import_indices(), ExcInternalError());
  ArrayView<unsigned int, MemorySpace::Default> temp_array_view_x(
    temp_array_view.data(), x.n_import_indices());
  {
    Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghosts(
      "ghosts", ghost_array_view.size());
    ArrayView<unsigned int, MemorySpace::Default> ghosts_view(ghosts.data(),
                                                              ghosts.size());
    Kokkos::deep_copy(ghosts, ghost_array);

    x.import_from_ghosted_array_start<unsigned int, MemorySpace::Default>(
      VectorOperation::add, 3, ghosts_view, temp_array_view_x, requests);
    x.import_from_ghosted_array_finish<unsigned int, MemorySpace::Default>(
      VectorOperation::add,
      temp_array_view_x,
      locally_owned_array_view,
      ghosts_view,
      requests);

    // check that the ghost entries are zeroed out in these calls
    deallog << "x ghost entries (should be zero up to index "
            << x.n_ghost_indices() << "):" << std::endl;
    print_device_view(ghosts);
  }
  deallog << "From reduced ghosts 2: ";
  print_device_view(locally_owned_array);

  // now send a tight array from x and add into the existing entries
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghosts(
    "ghosts", x.n_ghost_indices());
  ArrayView<unsigned int, MemorySpace::Default> ghosts_view(ghosts.data(),
                                                            ghosts.size());
  Kokkos::deep_copy(ghosts, 1);

  x.import_from_ghosted_array_start<unsigned int, MemorySpace::Default>(
    VectorOperation::add, 3, ghosts_view, temp_array_view_x, requests);
  x.import_from_ghosted_array_finish<unsigned int, MemorySpace::Default>(
    VectorOperation::add,
    temp_array_view_x,
    locally_owned_array_view,
    ghosts_view,
    requests);
  deallog << "From tight reduced ghosts 2: ";
  print_device_view(locally_owned_array);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;
  test();
}
