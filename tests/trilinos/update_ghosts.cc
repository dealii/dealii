// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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
