// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check SparsityTools::distribute_sparsity_pattern for BlockDynamicSP

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"



void
test_mpi()
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());

  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  unsigned int num_local = 10;
  unsigned int n         = numprocs * num_local;
  IndexSet     locally_owned_dofs(n);
  locally_owned_dofs.add_range(myid * num_local, (myid + 1) * num_local);

  IndexSet locally_rel(locally_owned_dofs);
  if (myid > 0)
    locally_rel.add_range((myid - 1) * num_local, (myid + 0) * num_local);
  if (myid < numprocs - 1)
    locally_rel.add_range((myid + 1) * num_local, (myid + 2) * num_local);

  std::vector<IndexSet> partitioning;
  partitioning.push_back(locally_rel);

  BlockDynamicSparsityPattern csp(partitioning);

  for (unsigned int i = 0; i < n; ++i)
    csp.add(i, myid);

  if (myid == 0)
    {
      deallog << "blocks: " << csp.n_block_rows() << 'x' << csp.n_block_cols()
              << std::endl;
      deallog << "size: " << csp.n_rows() << 'x' << csp.n_cols() << std::endl;
    }

  SparsityTools::distribute_sparsity_pattern(csp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_rel);
  /*  {
      std::ofstream
     f((std::string("after")+Utilities::int_to_string(myid)).c_str());
      csp.print(f);
      }*/

  // checking...
  for (unsigned int r = 0; r < num_local; ++r)
    {
      unsigned int indx = r + myid * num_local;
      unsigned int len  = csp.row_length(indx);

      // std::cout << "myid=" << myid << " idx=" << indx << " len=" << len
      // <<std::endl;

      if (myid > 0 && myid < numprocs - 1)
        Assert(len == 3, ExcInternalError());
      if (myid == 0 || myid == numprocs - 1)
        Assert(len == 2, ExcInternalError());

      Assert(csp.exists(indx, myid), ExcInternalError());
      if (myid > 0)
        Assert(csp.exists(indx, myid - 1), ExcInternalError());
      if (myid < numprocs - 1)
        Assert(csp.exists(indx, myid + 1), ExcInternalError());
    }


  // now a 2x2 block system where the 2,2 block has size 1x1:
  if (myid == 0)
    deallog << "part 2" << std::endl;

  IndexSet bla(1);
  if (myid == 0)
    bla.add_index(0);

  partitioning.push_back(bla);

  csp.reinit(partitioning);
  for (unsigned int i = 0; i < n; ++i)
    csp.add(i, myid);

  SparsityTools::distribute_sparsity_pattern(csp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_rel);

  if (myid == 0)
    {
      deallog << "blocks: " << csp.n_block_rows() << 'x' << csp.n_block_cols()
              << std::endl;
      deallog << "size: " << csp.n_rows() << 'x' << csp.n_cols() << std::endl;
    }

  // checking...
  for (unsigned int r = 0; r < num_local; ++r)
    {
      unsigned int indx = r + myid * num_local;
      unsigned int len  = csp.row_length(indx);

      // std::cout << "myid=" << myid << " idx=" << indx << " len=" << len
      // <<std::endl;

      if (myid > 0 && myid < numprocs - 1)
        Assert(len == 3, ExcInternalError());
      if (myid == numprocs - 1 || myid == 0)
        Assert(len == 2, ExcInternalError());

      Assert(csp.exists(indx, myid), ExcInternalError());
      if (myid > 0)
        Assert(csp.exists(indx, myid - 1), ExcInternalError());
      if (myid < numprocs - 1)
        Assert(csp.exists(indx, myid + 1), ExcInternalError());
    }

  if (myid == 0)
    deallog << "done" << std::endl;
}


int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test_mpi();
      deallog.pop();
    }
  else
    test_mpi();
}
