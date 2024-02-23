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



// allow operator= for BlockVectors if dest is same size or empty

#include <deal.II/base/index_set.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/generic_linear_algebra.h>

#include <iostream>
#include <vector>

#include "../tests.h"

#include "gla.h"

template <class LA>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  IndexSet block1(10);
  if (myid == 0)
    block1.add_range(0, 7);
  if (myid == 1)
    block1.add_range(7, 10);

  IndexSet block2(numproc);
  block2.add_index(myid);

  std::vector<IndexSet> partitioning;
  partitioning.push_back(block1);
  typename LA::MPI::BlockVector v_1(partitioning, MPI_COMM_WORLD);
  partitioning.push_back(block2);

  std::vector<IndexSet> relevant = partitioning;
  relevant[0].add_index(0);
  relevant[1].add_range(0, numproc);

  typename LA::MPI::BlockVector v_2(partitioning, MPI_COMM_WORLD);

  {
    typename LA::MPI::BlockVector x(partitioning, MPI_COMM_WORLD);
    x = v_2;
  }
  {
    typename LA::MPI::BlockVector x;
    x = v_2;
  }
  {
    deal_II_exceptions::disable_abort_on_exception();
    try
      {
        typename LA::MPI::BlockVector x = v_1;
        x                               = v_2; // error
      }
    catch (const ExceptionBase &e)
      {
        deallog << "Exception: " << e.get_exc_name() << std::endl;
      }
  }

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  {
    deallog.push("PETSc");
    test<LA_PETSc>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos>();
    deallog.pop();
  }
}
