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



// block vectors

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
  if (numproc == 1)
    block1.add_range(0, 10);

  if (myid == 0)
    block1.add_range(0, 7);
  if (myid == 1)
    block1.add_range(7, 10);

  IndexSet block2(6);
  if (numproc == 1)
    block2.add_range(0, 6);

  if (myid == 0)
    block2.add_range(0, 2);
  if (myid == 1)
    block2.add_range(2, 6);

  std::vector<IndexSet> partitioning;
  partitioning.push_back(block1);
  partitioning.push_back(block2);

  typename LA::MPI::BlockVector v(partitioning, MPI_COMM_WORLD);
  v = 0.1;

  deallog << "l2norm: " << v.l2_norm() << std::endl;

  v(myid) = myid;
  v.compress(VectorOperation::insert);

  deallog << "size: " << v.size() << std::endl;
  deallog << "block(0).size: " << v.block(0).size() << std::endl;
  deallog << "block(1).size: " << v.block(1).size() << std::endl;
  if (block1.n_elements() > 0)
    deallog << "my first entry: "
            << get_real_assert_zero_imag(v(block1.nth_index_in_set(0)))
            << std::endl;

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
