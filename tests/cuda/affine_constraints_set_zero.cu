// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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



// check AffineConstraints<double>::set_zero(Vector) for parallel distributed
// vectors

#include <deal.II/base/cuda_size.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"


__global__ void
initialize_vector(double *vector, int local_size, int offset)
{
  const int index = threadIdx.x + blockIdx.x * blockDim.x;
  if (index < local_size)
    vector[index] = 1.0 + index + offset;
}


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  IndexSet local_active;
  local_active.set_size(2 * numproc);
  local_active.add_range(myid * numproc, (myid + 1) * numproc);

  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> v;
  v.reinit(local_active, complete_index_set(2 * numproc), MPI_COMM_WORLD);

  const int n_blocks = 1 + v.size() / CUDAWrappers::block_size;
  initialize_vector<<<n_blocks, CUDAWrappers::block_size>>>(v.get_values(),
                                                            numproc,
                                                            myid * numproc);
  v.compress(VectorOperation::insert);

  AffineConstraints<double> cm;
  cm.add_line(numproc * myid + 1);
  cm.close();

  deallog << "vector before:" << std::endl;
  v.print(deallog.get_file_stream());

  deallog << std::endl;
  deallog << "CM:" << std::endl;
  cm.print(deallog.get_file_stream());

  cm.set_zero(v);

  deallog << "vector after:" << std::endl;
  v.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  init_cuda();

  test();
  return 0;
}
