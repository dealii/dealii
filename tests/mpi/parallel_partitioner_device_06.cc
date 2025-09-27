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
// regarding the export_to_ghosted_array() calls
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

  const unsigned int set = 200;
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
  local_relevant_1 = local_owned;

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

  // set up a locally owned array with some entries
  std::vector<unsigned int> cpu_locally_owned_data(local_size);
  for (unsigned int i = 0; i < local_size; ++i)
    cpu_locally_owned_data[i] = my_start + i;
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space>
    locally_owned_data("locally_owned_data", local_size);
  ArrayView<unsigned int, MemorySpace::Default> locally_owned_data_view(
    locally_owned_data.data(), local_size);
  Kokkos::deep_copy(
    locally_owned_data,
    Kokkos::View<unsigned *, Kokkos::HostSpace>(cpu_locally_owned_data.data(),
                                                cpu_locally_owned_data.size()));

  // set up a ghost array
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghosts(
    "ghosts", v.n_ghost_indices());
  ArrayView<unsigned int, MemorySpace::Default> ghosts_view(ghosts.data(),
                                                            ghosts.size());
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> temp_array(
    "temp_array", v.n_import_indices());
  ArrayView<unsigned int, MemorySpace::Default> temp_array_view(
    temp_array.data(), temp_array.size());

  std::vector<MPI_Request> requests;

  // send the full array
  v.export_to_ghosted_array_start<unsigned int, MemorySpace::Default>(
    3, locally_owned_data_view, temp_array_view, ghosts_view, requests);
  v.export_to_ghosted_array_finish<unsigned int, MemorySpace::Default>(
    ghosts_view, requests);
  deallog << "All ghosts: ";
  print_device_view(ghosts);

  // send only the array in w
  Kokkos::deep_copy(ghosts, 0);

  Assert(temp_array_view.size() >= w.n_import_indices(), ExcInternalError());
  ArrayView<unsigned int, MemorySpace::Default> temp_array_view_w(
    temp_array_view.data(), w.n_import_indices());
  w.export_to_ghosted_array_start<unsigned int, MemorySpace::Default>(
    3, locally_owned_data_view, temp_array_view_w, ghosts_view, requests);

  // start a second send operation for the x partitioner in parallel to make
  // sure communication does not get messed up
  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> temp_array2(
    "temp_array2", x.n_import_indices());
  ArrayView<unsigned int, MemorySpace::Default> temp_array2_view(
    temp_array2.data(), temp_array2.size());

  Kokkos::View<unsigned *, MemorySpace::Default::kokkos_space> ghosts2(
    "ghosts2", x.n_ghost_indices());
  ArrayView<unsigned int, MemorySpace::Default> ghosts2_view(ghosts2.data(),
                                                             ghosts2.size());

  std::vector<MPI_Request> requests2;
  x.export_to_ghosted_array_start<unsigned int, MemorySpace::Default>(
    4, locally_owned_data_view, temp_array2_view, ghosts2_view, requests2);

  w.export_to_ghosted_array_finish<unsigned int, MemorySpace::Default>(
    ghosts_view, requests);
  deallog << "Ghosts on reduced 1: ";
  print_device_view(ghosts);

  Kokkos::deep_copy(ghosts, 0);

  Assert(temp_array_view.size() >= x.n_import_indices(), ExcInternalError());
  ArrayView<unsigned int, MemorySpace::Default> temp_array_view_x(
    temp_array_view.data(), x.n_import_indices());
  x.export_to_ghosted_array_start<unsigned int, MemorySpace::Default>(
    3, locally_owned_data_view, temp_array_view_x, ghosts_view, requests);
  x.export_to_ghosted_array_finish<unsigned int, MemorySpace::Default>(
    ghosts_view, requests);
  deallog << "Ghosts on reduced 2: ";
  print_device_view(ghosts);

  x.export_to_ghosted_array_finish<unsigned int, MemorySpace::Default>(
    ghosts2_view, requests2);
  deallog << "Ghosts on reduced 2 without excess entries: ";
  print_device_view(ghosts2);

  x.export_to_ghosted_array_start<unsigned int, MemorySpace::Default>(
    3, locally_owned_data_view, temp_array_view_x, ghosts_view, requests);
  x.export_to_ghosted_array_finish<unsigned int, MemorySpace::Default>(
    ghosts_view, requests);
  deallog << "Ghosts on reduced 2: ";
  print_device_view(ghosts);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;
  test();
}
