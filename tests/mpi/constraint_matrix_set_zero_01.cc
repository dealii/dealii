// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double>::set_zero(Vector) for parallel vectors.
// this documents a bug introduced in r29678 for block vectors
// that was fixed in 29940.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"



void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  std::vector<IndexSet> local_active(2);

  // block 0:
  local_active[0].set_size(numproc);
  local_active[0].add_range(myid, myid + 1);

  // block 1:
  local_active[1].set_size(2 * numproc);
  local_active[1].add_range(myid * 2, myid * 2 + 2);

  PETScWrappers::MPI::BlockVector v(local_active, MPI_COMM_WORLD);

  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = 1.0 + i;
  v.compress(VectorOperation::insert);

  IndexSet local_active_together(3 * numproc);
  local_active_together.add_range(myid, myid + 1);
  local_active_together.add_range(numproc + myid * 2, numproc + myid * 2 + 2);

  AffineConstraints<PetscScalar> cm(local_active_together,
                                    local_active_together);
  cm.constrain_dof_to_zero(numproc + myid * 2);
  cm.close();

  deallog << "vector before:" << std::endl;
  v.print(deallog.get_file_stream());

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

  test();
  return 0;
}
