// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
