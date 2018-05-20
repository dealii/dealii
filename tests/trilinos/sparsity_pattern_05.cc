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

// Test copy_from(T)

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
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

  // Check copy_from(TrilinosWrappers::SparsityPattern):
  TrilinosWrappers::SparsityPattern B;
  B.copy_from(A);
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
  deallog << "OK" << std::endl;

  // copy_from(DynamicSparsityPattern)
  DynamicSparsityPattern dsp(4, 4);
  dsp.add(2, 3);
  B.copy_from(dsp);
  B.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;

  SparsityPattern sp(4, 4);
  sp.add(1, 2);
  sp.add(3, 3);
  sp.compress();
  B.copy_from(sp);
  B.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}
