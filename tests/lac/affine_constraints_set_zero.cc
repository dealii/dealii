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


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  IndexSet local_active;
  local_active.set_size(2 * numproc);
  local_active.add_range(myid * numproc, (myid + 1) * numproc);

  AffineConstraints<double> cm;
  cm.add_line(1);
  cm.add_line(2);
  cm.close();

  deallog << "CM:" << std::endl;
  cm.print(deallog.get_file_stream());

  using ExecutionSpace = MemorySpace::Default::kokkos_space::execution_space;
  ExecutionSpace exec;

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> ghosted;
  {
    ghosted.reinit(local_active,
                   complete_index_set(2 * numproc),
                   MPI_COMM_WORLD);
    auto ghosted_values = ghosted.get_values();

    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(exec, 0, numproc),
      KOKKOS_LAMBDA(int i) {
        int offset        = myid * numproc;
        ghosted_values[i] = 1.0 + i + offset;
      });
    ghosted.compress(VectorOperation::insert);

    deallog << "ghosted vector before:" << std::endl;
    ghosted.print(deallog.get_file_stream());

    cm.set_zero(ghosted);

    deallog << "ghosted vector after:" << std::endl;
    ghosted.print(deallog.get_file_stream());
  }

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> distributed;
  {
    distributed.reinit(local_active,
                       complete_index_set(2 * numproc),
                       MPI_COMM_WORLD);

    auto distributed_values = distributed.get_values();

    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(exec, 0, numproc),
      KOKKOS_LAMBDA(int i) {
        int offset            = myid * numproc;
        distributed_values[i] = 1.0 + i + offset;
      });
    exec.fence();

    distributed.compress(VectorOperation::insert);

    deallog << "distributed vector before:" << std::endl;
    distributed.print(deallog.get_file_stream());

    cm.set_zero(distributed);

    deallog << "distributed vector after:" << std::endl;
    distributed.print(deallog.get_file_stream());
  }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test();
  return 0;
}
