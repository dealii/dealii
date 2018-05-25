// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
#include <deal.II/lac/petsc_parallel_block_vector.h>

#include "../tests.h"

// Check that the base class vector and also the block vector class support
// update_ghost_values. This method doesn't do anything but is needed for
// genericity.

int
main(int argc, char **argv)
{
  using namespace dealii;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  PETScWrappers::VectorBase v;
  v.update_ghost_values();
  PETScWrappers::MPI::BlockVector bv;
  bv.update_ghost_values();
  deallog << "OK" << std::endl;
}
