// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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
