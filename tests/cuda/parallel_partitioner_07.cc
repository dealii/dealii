// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
print_cuda_view(const ArrayView<Number, MemorySpace::CUDA> cuda_view)
{
  std::vector<Number> cpu_values(cuda_view.size());
  Utilities::CUDA::copy_to_host(cuda_view.data(), cpu_values);
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
  std::vector<unsigned int> cpu_ghost_array(v.n_ghost_indices(), 1);
  std::unique_ptr<unsigned int[], void (*)(unsigned int *)> ghost_array(
    Utilities::CUDA::allocate_device_data<unsigned int>(cpu_ghost_array.size()),
    Utilities::CUDA::delete_device_data<unsigned int>);
  ArrayView<unsigned int, MemorySpace::CUDA> ghost_array_view(
    ghost_array.get(), cpu_ghost_array.size());
  Utilities::CUDA::copy_to_dev(cpu_ghost_array, ghost_array.get());

  // set up other arrays
  std::unique_ptr<unsigned int[], void (*)(unsigned int *)> locally_owned_array(
    Utilities::CUDA::allocate_device_data<unsigned int>(local_size),
    Utilities::CUDA::delete_device_data<unsigned int>);
  ArrayView<unsigned int, MemorySpace::CUDA> locally_owned_array_view(
    locally_owned_array.get(), local_size);

  std::unique_ptr<unsigned int[], void (*)(unsigned int *)> temp_array(
    Utilities::CUDA::allocate_device_data<unsigned int>(v.n_import_indices()),
    Utilities::CUDA::delete_device_data<unsigned int>);
  ArrayView<unsigned int, MemorySpace::CUDA> temp_array_view(
    temp_array.get(), v.n_import_indices());

  std::vector<MPI_Request> requests;

  // send the full array
  {
    std::unique_ptr<unsigned int[], void (*)(unsigned int *)> ghosts(
      Utilities::CUDA::allocate_device_data<unsigned int>(
        ghost_array_view.size()),
      Utilities::CUDA::delete_device_data<unsigned int>);
    ArrayView<unsigned int, MemorySpace::CUDA> ghosts_view(
      ghosts.get(), ghost_array_view.size());
    const cudaError_t cuda_error =
      cudaMemcpy(ghosts.get(),
                 ghost_array_view.data(),
                 ghost_array_view.size() * sizeof(unsigned int),
                 cudaMemcpyDeviceToDevice);
    AssertCuda(cuda_error);

    v.import_from_ghosted_array_start<unsigned int, MemorySpace::CUDA>(
      VectorOperation::add, 3, ghosts_view, temp_array_view, requests);
    v.import_from_ghosted_array_finish<unsigned int, MemorySpace::CUDA>(
      VectorOperation::add,
      temp_array_view,
      locally_owned_array_view,
      ghosts_view,
      requests);
    // check that the ghost entries are zeroed out in these calls
    deallog << "v ghost entries (should be zero up to index "
            << v.n_ghost_indices() - 1 << "):" << std::endl;
    print_cuda_view(ghosts_view);
  }
  deallog << "From all ghosts: ";
  print_cuda_view(locally_owned_array_view);

  // send only the array in w
  cudaError_t cuda_error =
    cudaMemset(locally_owned_array_view.data(),
               0,
               locally_owned_array_view.size() * sizeof(unsigned int));
  AssertCuda(cuda_error);
  Assert(temp_array_view.size() >= w.n_import_indices(), ExcInternalError());
  ArrayView<unsigned int, MemorySpace::CUDA> temp_array_view_w(
    temp_array_view.data(), w.n_import_indices());
  {
    std::unique_ptr<unsigned int[], void (*)(unsigned int *)> ghosts(
      Utilities::CUDA::allocate_device_data<unsigned int>(
        ghost_array_view.size()),
      Utilities::CUDA::delete_device_data<unsigned int>);
    ArrayView<unsigned int, MemorySpace::CUDA> ghosts_view(
      ghosts.get(), ghost_array_view.size());
    const cudaError_t cuda_error =
      cudaMemcpy(ghosts.get(),
                 ghost_array_view.data(),
                 ghost_array_view.size() * sizeof(unsigned int),
                 cudaMemcpyDeviceToDevice);
    AssertCuda(cuda_error);

    w.import_from_ghosted_array_start<unsigned int, MemorySpace::CUDA>(
      VectorOperation::add, 3, ghosts_view, temp_array_view_w, requests);
    w.import_from_ghosted_array_finish<unsigned int, MemorySpace::CUDA>(
      VectorOperation::add,
      temp_array_view_w,
      locally_owned_array_view,
      ghosts_view,
      requests);

    // check that the ghost entries are zeroed out in these calls
    deallog << "w ghost entries (should be zero up to index "
            << w.n_ghost_indices() - 1 << "):" << std::endl;
    print_cuda_view(ghosts_view);
  }
  deallog << "From reduced ghosts 1: ";
  print_cuda_view(locally_owned_array_view);

  // send only the array in x
  cuda_error =
    cudaMemset(locally_owned_array_view.data(),
               0,
               locally_owned_array_view.size() * sizeof(unsigned int));
  AssertCuda(cuda_error);
  Assert(temp_array_view.size() >= x.n_import_indices(), ExcInternalError());
  ArrayView<unsigned int, MemorySpace::CUDA> temp_array_view_x(
    temp_array_view.data(), x.n_import_indices());
  {
    std::unique_ptr<unsigned int[], void (*)(unsigned int *)> ghosts(
      Utilities::CUDA::allocate_device_data<unsigned int>(
        ghost_array_view.size()),
      Utilities::CUDA::delete_device_data<unsigned int>);
    ArrayView<unsigned int, MemorySpace::CUDA> ghosts_view(
      ghosts.get(), ghost_array_view.size());
    const cudaError_t cuda_error =
      cudaMemcpy(ghosts.get(),
                 ghost_array_view.data(),
                 ghost_array_view.size() * sizeof(unsigned int),
                 cudaMemcpyDeviceToDevice);
    AssertCuda(cuda_error);

    x.import_from_ghosted_array_start<unsigned int, MemorySpace::CUDA>(
      VectorOperation::add, 3, ghosts_view, temp_array_view_x, requests);
    x.import_from_ghosted_array_finish<unsigned int, MemorySpace::CUDA>(
      VectorOperation::add,
      temp_array_view_x,
      locally_owned_array_view,
      ghosts_view,
      requests);

    // check that the ghost entries are zeroed out in these calls
    deallog << "x ghost entries (should be zero up to index "
            << x.n_ghost_indices() << "):" << std::endl;
    print_cuda_view(ghosts_view);
  }
  deallog << "From reduced ghosts 2: ";
  print_cuda_view(locally_owned_array_view);

  // now send a tight array from x and add into the existing entries
  std::vector<unsigned int> cpu_ghosts(x.n_ghost_indices(), 1);
  std::unique_ptr<unsigned int[], void (*)(unsigned int *)> ghosts(
    Utilities::CUDA::allocate_device_data<unsigned int>(cpu_ghosts.size()),
    Utilities::CUDA::delete_device_data<unsigned int>);
  ArrayView<unsigned int, MemorySpace::CUDA> ghosts_view(ghosts.get(),
                                                         cpu_ghosts.size());
  Utilities::CUDA::copy_to_dev(cpu_ghosts, ghosts.get());

  x.import_from_ghosted_array_start<unsigned int, MemorySpace::CUDA>(
    VectorOperation::add, 3, ghosts_view, temp_array_view_x, requests);
  x.import_from_ghosted_array_finish<unsigned int, MemorySpace::CUDA>(
    VectorOperation::add,
    temp_array_view_x,
    locally_owned_array_view,
    ghosts_view,
    requests);
  deallog << "From tight reduced ghosts 2: ";
  print_cuda_view(locally_owned_array_view);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;
  init_cuda(true);
  test();
}
