// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// check Utilities::MPI::logical_and()

#include <deal.II/base/utilities.h>

#include "../tests.h"

void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  bool collective_false = Utilities::MPI::logical_and(false, MPI_COMM_WORLD);
  bool collective_alternating =
    Utilities::MPI::logical_and(((myid % 2) == 0) ? true : false,
                                MPI_COMM_WORLD);
  bool collective_true = Utilities::MPI::logical_and(true, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << std::boolalpha << collective_false << ' '
            << collective_alternating << ' ' << collective_true << std::endl;
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
      test();
      deallog.pop();
    }
  else
    test();
}
