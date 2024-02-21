// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
