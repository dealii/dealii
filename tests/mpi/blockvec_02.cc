// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test constructor/reinit of BlockVector with IndexSets and the conversion
// to a Vector

#include <deal.II/base/index_set.h>

#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

#include "../distributed_grids/coarse_grid_common.h"

void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  std::vector<IndexSet> local_active(2);

  // block 0:
  local_active[0].set_size(numproc);
  local_active[0].add_range(myid, myid + 1);

  local_active[1].set_size(2 * numproc);
  local_active[1].add_range(myid * 2, myid * 2 + 2);

  TrilinosWrappers::MPI::BlockVector v_block(local_active, MPI_COMM_WORLD);

  v_block(myid) = 100.0 + myid;

  v_block.block(1)(myid * 2)     = myid * 2.0;
  v_block.block(1)(myid * 2 + 1) = myid * 2.0 + 1.0;

  v_block.compress(VectorOperation::insert);

  deallog << "block size: " << v_block.size() << std::endl;
  deallog << "block size[0]: " << v_block.block(0).size() << std::endl;
  deallog << "block size[1]: " << v_block.block(1).size() << std::endl;
  v_block.locally_owned_elements().print(deallog.get_file_stream());

  TrilinosWrappers::MPI::Vector v;
  v.reinit(v_block);

  deallog << "size: " << v.size() << std::endl;
  v.locally_owned_elements().print(deallog.get_file_stream());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  test();
}
