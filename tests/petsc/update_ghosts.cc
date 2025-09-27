// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#include <deal.II/lac/petsc_block_vector.h>

#include "../tests.h"

// Check that the base class vector and also the block vector class support
// update_ghost_values. This method doesn't do anything but is needed for
// genericity.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  PETScWrappers::VectorBase v;
  v.update_ghost_values();
  PETScWrappers::MPI::BlockVector bv;
  bv.update_ghost_values();
  deallog << "OK" << std::endl;
}
