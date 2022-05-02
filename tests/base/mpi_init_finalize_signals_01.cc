// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Test signals in MPI_InitFinalize

#include <deal.II/base/mpi.h>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  initlog();

  Utilities::MPI::MPI_InitFinalize::signals.at_mpi_init.connect([=]() {
    deallog << "called via MPI_InitFinalize::signals.at_mpi_init" << std::endl;
  });

  Utilities::MPI::MPI_InitFinalize::signals.at_mpi_finalize.connect([=]() {
    deallog << "called via MPI_InitFinalize::signals.at_mpi_finalize"
            << std::endl;
  });

  deallog << "before MPI initialization" << std::endl;
  {
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    deallog << "after MPI initialization" << std::endl;
  }
  deallog << "after MPI finalization" << std::endl;
}
