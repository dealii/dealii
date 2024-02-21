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

// test MPI::Partitioner update_ghosts() and compress() in case we have
// empty owned DoFs

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>

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

template <typename Number = double>
void
test()
{
  const unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // setup index sets
  //                            subset:                    is2
  //                            ghost:       8 9  10 11    is3
  //      rank 0 :  00 01 02 03 04 05 06 07 00 00 00 00
  //      rank 1 :  00 00 00 00 00 00 00 00 00 00 00 00
  //                            ghost:      0  1  2  3     is3
  //                            subset:        1  2        is2
  //
  // expected result update ghosts()
  //
  //      rank 0 :  00 01 02 03 04 05 06 07 00 00 00 00
  //      rank 1 :  00 00 00 00 00 00 00 00 00 01 02 00
  //
  // compress(insert) -- does not change anything but zero ghosts
  //
  // set rank1 ghosts to: 00 10 20 00
  // compress(add)
  //
  //      rank 0 :  00 11 22 03 04 05 06 07 00 00 00 00
  //      rank 1 :  00 00 00 00 00 00 00 00 00 10 20 00


  IndexSet is1(16), is2(16), is3(16);

  if (rank == 0)
    {
      is1.add_range(0, 8);
      // note: empty is2
      is3.add_range(8, 12);
    }
  else if (rank == 1)
    {
      is1.add_range(8, 16);
      is2.add_index(1);
      is2.add_index(2);
      is3.add_range(0, 4);
    }

  // create partitioner
  std::shared_ptr<Utilities::MPI::Partitioner> partitioner(
    new Utilities::MPI::Partitioner(is1, MPI_COMM_WORLD));
  partitioner->set_ghost_indices(is3);
  std::shared_ptr<Utilities::MPI::Partitioner> tight_partitioner(
    new Utilities::MPI::Partitioner(is1, MPI_COMM_WORLD));
  tight_partitioner->set_ghost_indices(is2, is3);

  // create vector
  std::vector<Number> cpu_owned(rank == 0 ? 8 : 0);
  for (unsigned int i = 0; i < cpu_owned.size(); ++i)
    cpu_owned[i] = i;
  Kokkos::View<Number *, MemorySpace::Default::kokkos_space> owned(
    "owned", rank == 0 ? 8 : 0);
  ArrayView<Number, MemorySpace::Default>             owned_view(owned.data(),
                                                     owned.size());
  MemorySpace::Default::kokkos_space::execution_space exec;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<decltype(exec)>(exec, 0, owned.size()),
    KOKKOS_LAMBDA(int i) { owned(i) = i; });
  exec.fence();

  Kokkos::View<Number *, MemorySpace::Default::kokkos_space> ghost("ghost", 4);
  ArrayView<Number, MemorySpace::Default> ghost_view(ghost.data(),
                                                     ghost.size());
  Kokkos::deep_copy(ghost, 0);

  // update ghost values
  // vector of requests
  std::vector<MPI_Request> requests;
  std::vector<MPI_Request> compress_requests;

  // allocate temporal array
  Kokkos::View<Number *, MemorySpace::Default::kokkos_space> tmp_data(
    "tmp_data", tight_partitioner->n_import_indices());
  ArrayView<Number, MemorySpace::Default> tmp_data_view(tmp_data.data(),
                                                        tmp_data.size());

  // begin exchange, and ...
  tight_partitioner
    ->export_to_ghosted_array_start<Number, MemorySpace::Default>(
      0, owned_view, tmp_data_view, ghost_view, requests);

  // ... finish exchange
  tight_partitioner
    ->export_to_ghosted_array_finish<Number, MemorySpace::Default>(ghost_view,
                                                                   requests);

  auto print = [&]() {
    deallog << "owned:" << std::endl;
    print_device_view(owned);
    deallog << "ghost:" << std::endl;
    print_device_view(ghost);
  };

  deallog << "update ghosts()" << std::endl;
  print();

  Kokkos::View<Number *, MemorySpace::Default::kokkos_space> import_data(
    "import_data", tight_partitioner->n_import_indices());
  ArrayView<Number, MemorySpace::Default> import_data_view(tmp_data.data(),
                                                           import_data.size());

  // now do insert:
  auto compress = [&](VectorOperation::values operation) {
    const unsigned int counter = 0;
    tight_partitioner
      ->import_from_ghosted_array_start<Number, MemorySpace::Default>(
        operation, counter, ghost_view, import_data_view, compress_requests);

    tight_partitioner
      ->import_from_ghosted_array_finish<Number, MemorySpace::Default>(
        operation, import_data_view, owned_view, ghost_view, compress_requests);
  };

  deallog << "compress(insert)" << std::endl;
  compress(VectorOperation::insert);
  print();

  if (rank == 1)
    {
      Kokkos::deep_copy(Kokkos::subview(ghost, 1), 10);
      Kokkos::deep_copy(Kokkos::subview(ghost, 2), 20);
    }

  deallog << "compress(add)" << std::endl;
  compress(VectorOperation::add);
  print();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();

  return 0;
}
