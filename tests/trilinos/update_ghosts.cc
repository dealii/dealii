// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include "../tests.h"


// Check that the block vector class support update_ghost_values. This method
// doesn't do anything but is needed for genericity.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  TrilinosWrappers::MPI::BlockVector bv;
  bv.update_ghost_values();
  deallog << "OK" << std::endl;
}
