// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

template <typename Number>
void
print_cuda_view(const ArrayView<Number, MemorySpace::CUDA> cuda_view)
{
  std::vector<Number> cpu_values(cuda_view.size());
  Utilities::CUDA::copy_to_host(cuda_view.data(), cpu_values);
  for (Number value : cpu_values)
    deallog << value << " ";
  deallog << std::endl;
}

__global__ void
set_value(double *values_dev, unsigned int index, double val)
{
  values_dev[index] = val;
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
  std::unique_ptr<Number[], void (*)(Number *)> owned(
    Utilities::CUDA::allocate_device_data<Number>(cpu_owned.size()),
    Utilities::CUDA::delete_device_data<Number>);
  ArrayView<Number, MemorySpace::CUDA> owned_view(owned.get(),
                                                  cpu_owned.size());
  Utilities::CUDA::copy_to_dev(cpu_owned, owned.get());

  std::vector<Number>                           cpu_ghost(4, 0);
  std::unique_ptr<Number[], void (*)(Number *)> ghost(
    Utilities::CUDA::allocate_device_data<Number>(cpu_ghost.size()),
    Utilities::CUDA::delete_device_data<Number>);
  ArrayView<Number, MemorySpace::CUDA> ghost_view(ghost.get(),
                                                  cpu_ghost.size());
  Utilities::CUDA::copy_to_dev(cpu_ghost, ghost.get());

  // update ghost values
  // vector of requests
  std::vector<MPI_Request> requests;
  std::vector<MPI_Request> compress_requests;

  // allocate temporal array
  std::unique_ptr<Number[], void (*)(Number *)> tmp_data(
    Utilities::CUDA::allocate_device_data<Number>(
      tight_partitioner->n_import_indices()),
    Utilities::CUDA::delete_device_data<Number>);
  ArrayView<Number, MemorySpace::CUDA> tmp_data_view(
    tmp_data.get(), tight_partitioner->n_import_indices());

  // begin exchange, and ...
  tight_partitioner->export_to_ghosted_array_start<Number, MemorySpace::CUDA>(
    0, owned_view, tmp_data_view, ghost_view, requests);

  // ... finish exchange
  tight_partitioner->export_to_ghosted_array_finish<Number, MemorySpace::CUDA>(
    ghost_view, requests);

  auto print = [&]() {
    deallog << "owned:" << std::endl;
    print_cuda_view(owned_view);
    deallog << "ghost:" << std::endl;
    print_cuda_view(ghost_view);
  };

  deallog << "update ghosts()" << std::endl;
  print();

  std::unique_ptr<Number[], void (*)(Number *)> import_data(
    Utilities::CUDA::allocate_device_data<Number>(
      tight_partitioner->n_import_indices()),
    Utilities::CUDA::delete_device_data<Number>);
  ArrayView<Number, MemorySpace::CUDA> import_data_view(
    tmp_data.get(), tight_partitioner->n_import_indices());

  // now do insert:
  auto compress = [&](VectorOperation::values operation) {
    const unsigned int counter = 0;
    tight_partitioner
      ->import_from_ghosted_array_start<Number, MemorySpace::CUDA>(
        operation, counter, ghost_view, import_data_view, compress_requests);

    tight_partitioner
      ->import_from_ghosted_array_finish<Number, MemorySpace::CUDA>(
        operation, import_data_view, owned_view, ghost_view, compress_requests);
  };

  deallog << "compress(insert)" << std::endl;
  compress(VectorOperation::insert);
  print();

  if (rank == 1)
    {
      set_value<<<1, 1>>>(ghost.get(), 1, 10);
      set_value<<<1, 1>>>(ghost.get(), 2, 20);
    }

  deallog << "compress(add)" << std::endl;
  compress(VectorOperation::add);
  print();
}

int
main(int argc, char **argv)
{
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  init_cuda(true);

  test();

  return 0;
}
