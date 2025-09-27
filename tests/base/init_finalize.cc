// ---------------------------------------------------------------------
//
// Copyright (C) 2024 by the deal.II authors
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
