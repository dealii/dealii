// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Try to compile a file that includes PSBLAS headers, and do an MPI hello world
// program using psblas interfaces.

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/psblas_precondition.h>

#include "../tests.h"


using namespace dealii;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm                         mpi_communicator = MPI_COMM_WORLD;

  mpi_initlog();
  // create and free an amg preconditioner
  amg_c_dprec *psblas_preconditioner = amg_c_dprec_new();
  free(psblas_preconditioner);
  deallog << "Detection AMG4PSBLAS: OK" << std::endl;
}
