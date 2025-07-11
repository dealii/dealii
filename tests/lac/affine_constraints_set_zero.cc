// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double>::set_zero(Vector) for parallel distributed
// vectors

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
  cm.constrain_dof_to_zero(1);
  cm.constrain_dof_to_zero(2);
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
