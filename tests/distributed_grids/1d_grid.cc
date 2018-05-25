// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



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
