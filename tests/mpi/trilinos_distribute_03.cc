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



// document bug in the new AffineConstraints<double>::distribute() from r29593
// if one vector owns all DoFs, then distribute() hangs.

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"



void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // create a vector that consists of elements indexed from 0 to 100
  // everything belongs to proc 0
  IndexSet idx(100);
  if (myid == 0)
    idx.add_range(0, 100);

  TrilinosWrappers::MPI::Vector vec(idx, MPI_COMM_WORLD);

  for (unsigned int i = vec.local_range().first; i < vec.local_range().second;
       ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  AffineConstraints<double> cm;
  cm.close();

  cm.distribute(vec);

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
