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

// Test the move constructor of the sparsity pattern and check the sparsity
// pattern whose data have been stolen can be reinitialized

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

int
main(int argc, char** argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  IndexSet partitioning(3);

  partitioning.add_range(0, 3);

  // Add element (2,1) to the matrix
  TrilinosWrappers::SparsityPattern A(partitioning);
  A.add(2, 1);
  A.compress();

  // Check the move constructor
  TrilinosWrappers::SparsityPattern B(std::move(A));
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      {
        if((i == 2) && (j == 1))
          {
            AssertThrow(B.exists(i, j) == true, ExcInternalError());
          }
        else
          {
            AssertThrow(B.exists(i, j) == false, ExcInternalError());
          }
      }

  // Check that A can be reinitialized
  A.reinit(partitioning);
  A.add(1, 2);
  A.compress();
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      {
        if((i == 1) && (j == 2))
          {
            AssertThrow(A.exists(i, j) == true, ExcInternalError());
          }
        else
          {
            AssertThrow(A.exists(i, j) == false, ExcInternalError());
          }
      }

  deallog << "OK" << std::endl;
}
