// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check that InitFinalize can initialize/finalize only the selected library

#include <deal.II/base/init_finalize.h>

#include <Kokkos_Core.hpp>
#include <p4est_vtk.h>

#include "../tests.h"

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  mpi_initlog();
  {
    // Only initialize P4EST. If it tries to initialize MPI or Kokkos, we will
    // get an error.
    InitFinalize init_finalize(argc, argv, InitializeLibrary::P4EST);
    // p4est_is_initialized was introduced in version 2.8. Instead we create the
    // connectivity. If P4EST was not initialized the code below will crash.
    auto *conn = p4est_connectivity_new_unitsquare();
    p4est_connectivity_destroy(conn);
  }
  deallog << "OK" << std::endl;

  Kokkos::finalize();

  MPI_Finalize();

  return 0;
}
