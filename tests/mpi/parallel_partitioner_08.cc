// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

// test MPI::Partitioner update_ghosts() and compress() in case we have
// empty owned DoFs

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>

#include "../tests.h"

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
  AlignedVector<Number> owned(rank == 0 ? 8 : 0);
  AlignedVector<Number> ghost(4);

  for (unsigned int i = 0; i < 4; ++i)
    ghost[i] = 0.;

  if (rank == 0)
    for (int i = 0; i < 8; i++)
      owned[i] = i;

  // update ghost values
  // vector of requests
  std::vector<MPI_Request> requests;
  std::vector<MPI_Request> compress_requests;

  // allocate temporal array
  AlignedVector<Number> tmp_data;
  tmp_data.resize_fast(tight_partitioner->n_import_indices());

  // begin exchange, and ...
  tight_partitioner->export_to_ghosted_array_start(
    0,
    ArrayView<const Number>(owned.begin(), owned.size()),
    ArrayView<Number>(tmp_data.begin(), tight_partitioner->n_import_indices()),
    ArrayView<Number>(ghost.begin(), ghost.size()),
    requests);

  // ... finish exchange
  tight_partitioner->export_to_ghosted_array_finish(
    ArrayView<Number>(ghost.begin(), ghost.size()), requests);

  auto print = [&]() {
    deallog << "owned:" << std::endl;
    for (auto el : owned)
      deallog << el << " ";
    deallog << std::endl << "ghost:" << std::endl;
    for (auto el : ghost)
      deallog << el << " ";
    deallog << std::endl;
  };

  deallog << "update ghosts()" << std::endl;
  print();

  AlignedVector<Number> import_data;
  import_data.resize_fast(tight_partitioner->n_import_indices());

  // now do insert:
  auto compress = [&](VectorOperation::values operation) {
    const unsigned int counter = 0;
    tight_partitioner->import_from_ghosted_array_start(
      operation,
      counter,
      ArrayView<Number>(ghost.begin(), ghost.size()),
      ArrayView<Number>(import_data.begin(),
                        tight_partitioner->n_import_indices()),
      compress_requests);

    tight_partitioner->import_from_ghosted_array_finish(
      operation,
      ArrayView<const Number>(import_data.begin(),
                              tight_partitioner->n_import_indices()),
      ArrayView<Number>(owned.begin(), owned.size()),
      ArrayView<Number>(ghost.begin(), ghost.size()),
      compress_requests);
  };

  deallog << "compress(insert)" << std::endl;
  compress(VectorOperation::insert);
  print();

  if (rank == 1)
    {
      ghost[1] = 10;
      ghost[2] = 20;
    }

  deallog << "compress(add)" << std::endl;
  compress(VectorOperation::add);
  print();
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  test();

  return 0;
}
