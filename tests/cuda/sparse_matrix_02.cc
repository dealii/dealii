// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Check print and print_format

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/lac/cuda_sparse_matrix.h>

#include "../tests.h"

#include "../testmatrix.h"


void
test(Utilities::CUDA::Handle &cuda_handle)
{
  // Build the sparse matrix on the host
  const unsigned int size = 3;
  unsigned int       dim  = (size - 1) * (size - 1);

  FDMatrix        testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double> A(structure);
  testproblem.upwind(A, true);
  A.print(deallog.get_file_stream());
  A.print_formatted(deallog.get_file_stream());

  // Create the sparse matrix on the device
  CUDAWrappers::SparseMatrix<double> A_dev(cuda_handle, A);
  A_dev.print(deallog.get_file_stream());
  A_dev.print_formatted(deallog.get_file_stream());
}

int
main()
{
  initlog();
  deallog.depth_console(0);

  init_cuda();

  Utilities::CUDA::Handle cuda_handle;

  test(cuda_handle);

  deallog << "OK" << std::endl;

  return 0;
}
