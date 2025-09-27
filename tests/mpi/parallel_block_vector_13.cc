// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check ghost handling on parallel block vectors for large
// number of blocks with split compress_start()/update_ghosts_start().
// almost copy-paste of parallel_block_vector_02.cc

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;


  // each processor from processor 1 to 8 owns 2 indices (the other processors
  // do not own any dof), and all processors are ghosting element 1
  IndexSet local_owned(std::min(16U, numproc * 2));
  if (myid < 8)
    local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);

  LinearAlgebra::distributed::Vector<double> v(local_owned,
                                               local_relevant,
                                               MPI_COMM_WORLD);

  // set local values
  if (myid < 8)
    {
      v(myid * 2)     = myid * 2.0;
      v(myid * 2 + 1) = myid * 2.0 + 1.0;
    }
  v.compress(VectorOperation::insert);

  const unsigned int n_blocks = 107;

  LinearAlgebra::distributed::BlockVector<double> w(n_blocks);
  for (unsigned int i = 0; i < n_blocks; ++i)
    {
      w.block(i) = v;
      w.block(i) *= (i + 1);
    }
  w.collect_sizes();

  // see if we can access the ghosted entry 1 in each block with the correct
  // value: the initialization should non have updated the ghost values, so
  // all other processors except 0 should have zero entry on the global index
  // 1
  for (unsigned int i = 0; i < n_blocks; ++i)
    if (myid == 0)
      {
        AssertDimension(i + 1, (unsigned int)w.block(i)(1));
      }
    else
      {
        AssertDimension(0, (unsigned int)w.block(i)(1));
      }

  // import ghost values, all processors should still have i+1
  w.update_ghost_values();
  for (unsigned int i = 0; i < n_blocks; ++i)
    AssertDimension(i + 1, (unsigned int)w.block(i)(1));

  // zero out ghosts, now all processors except processor 1 should have 0.
  w.zero_out_ghost_values();
  for (unsigned int i = 0; i < n_blocks; ++i)
    if (myid == 0)
      {
        AssertDimension(i + 1, (unsigned int)w.block(i)(1));
      }
    else
      {
        AssertDimension(0, (unsigned int)w.block(i)(1));
      }

  // create a vector copy that gets the entries from w. First, it should not
  // have updated the ghosts because it is created from an empty state.
  LinearAlgebra::distributed::BlockVector<double> x(w);
  Assert(x.has_ghost_elements() == false, ExcInternalError());
  for (unsigned int i = 0; i < n_blocks; ++i)
    if (myid == 0)
      {
        AssertDimension(i + 1, (unsigned int)x.block(i)(1));
      }
    else
      {
        AssertDimension(0, (unsigned int)x.block(i)(1));
      }

  // now we zero the vector, which should disable ghost elements
  x = 0;
  Assert(x.has_ghost_elements() == false, ExcInternalError());

  // we copy from w (i.e., the same vector but one that does not have ghosts
  // enabled) -> should not have ghosts enabled
  x = w;
  Assert(x.has_ghost_elements() == false, ExcInternalError());
  for (unsigned int i = 0; i < n_blocks; ++i)
    if (myid == 0)
      {
        AssertDimension(i + 1, (unsigned int)x.block(i)(1));
      }
    else
      {
        AssertDimension(0, (unsigned int)x.block(i)(1));
      }

  x.update_ghost_values();
  Assert(x.has_ghost_elements() == true, ExcInternalError());


  // add something to entry 1 on all processors
  w(1) += myid + 1;
  w.compress(VectorOperation::add);
  if (myid == 0)
    AssertDimension((unsigned int)w(1), 1 + (numproc * (numproc + 1)) / 2);

  // add again and check if everything is still correct
  w(1 + v.size()) += myid + 1;
  w.compress(VectorOperation::add);
  if (myid == 0)
    AssertDimension((unsigned int)w(1), 1 + (numproc * (numproc + 1)) / 2);
  if (myid == 0)
    AssertDimension((unsigned int)w(v.size() + 1),
                    2 + (numproc * (numproc + 1)) / 2);

  w.update_ghost_values();
  AssertDimension((unsigned int)w(1), 1 + (numproc * (numproc + 1)) / 2);
  AssertDimension((unsigned int)w(v.size() + 1),
                  2 + (numproc * (numproc + 1)) / 2);

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
