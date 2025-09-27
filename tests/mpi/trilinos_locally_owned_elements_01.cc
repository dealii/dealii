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



// test TrilinosVector::locally_owned_elements


#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // create non-contiguous index set
  {
    AssertThrow(n_processes == 2, ExcNotImplemented());
    IndexSet index(10);
    for (unsigned int i = 0; i < 10; i += 2)
      index.add_range(i + myid, i + myid + 1);
    index.compress();

    TrilinosWrappers::MPI::Vector vec(index, MPI_COMM_WORLD);

    IndexSet index2 = vec.locally_owned_elements();
    AssertThrow(index == index2, ExcInternalError());
  }

  // create contiguous index set
  {
    IndexSet index(10);
    index.add_range(5 * myid, 5 * (myid + 1));
    index.compress();

    TrilinosWrappers::MPI::Vector vec(index, MPI_COMM_WORLD);

    IndexSet index2 = vec.locally_owned_elements();
    AssertThrow(index == index2, ExcInternalError());
  }

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

      test<2>();
    }
  else
    {
      test<2>();
    }
}
