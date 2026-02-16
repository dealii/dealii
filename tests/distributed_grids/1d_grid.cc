// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test that we can compile p::d::Tria<1>

#include <deal.II/distributed/tria.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  deal_II_exceptions::disable_abort_on_exception();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  try
    {
      parallel::distributed::Triangulation<1> tr(MPI_COMM_WORLD);
    }
  catch (const std::exception &exc)
    {
      deallog << "This test has to throw an exception" << std::endl;
    }
}
