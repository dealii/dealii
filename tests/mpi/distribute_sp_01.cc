// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check SparsityTools::distribute_sparsity_pattern with the slow part of
// supplying the number of rows on all processes

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

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

  unsigned int                         num_local = 10;
  unsigned int                         n         = numprocs * num_local;
  std::vector<types::global_dof_index> rows_per_cpu;
  for (unsigned int i = 0; i < numprocs; ++i)
    rows_per_cpu.push_back(num_local);

  IndexSet locally_rel(n);
  locally_rel.add_range(myid * num_local, (myid + 1) * num_local);
  if (myid > 0)
    locally_rel.add_range((myid - 1) * num_local, (myid + 0) * num_local);
  if (myid < numprocs - 1)
    locally_rel.add_range((myid + 1) * num_local, (myid + 2) * num_local);

  DynamicSparsityPattern csp(n, n, locally_rel);

  for (unsigned int i = 0; i < n; ++i)
    csp.add(i, myid);

  SparsityTools::distribute_sparsity_pattern(csp,
                                             rows_per_cpu,
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
