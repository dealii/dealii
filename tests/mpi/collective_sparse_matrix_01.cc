// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::MPI::sum() for sparse matrix

#include <deal.II/base/mpi.h>

#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"

#include "../testmatrix.h"



void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int size = 50;

  FDMatrix     testproblem(size, size);
  unsigned int dim = (size - 1) * (size - 1);

  SparsityPattern sparsity(dim, dim, size);
  testproblem.five_point_structure(sparsity);
  sparsity.compress();

  SparseMatrix<double> matrix(sparsity);
  const double         val = std::pow(10, myid);
  for (SparsityPattern::const_iterator it = sparsity.begin();
       it != sparsity.end();
       ++it)
    {
      const auto i = (*it).row();
      const auto j = (*it).column();
      matrix.add(i, j, -val);
      matrix.add(i, i, val);
    }

  // compare with FullMatrix:
  FullMatrix<double> full(dim, dim);
  full.copy_from(matrix);

  // deallog << "Local:" << std::endl;
  // matrix.print(deallog.get_file_stream());
  Utilities::MPI::sum(matrix, MPI_COMM_WORLD, matrix);
  // deallog << "Global:" << std::endl;
  // matrix.print(deallog.get_file_stream());

  Utilities::MPI::sum(full, MPI_COMM_WORLD, full);

  for (SparsityPattern::const_iterator it = sparsity.begin();
       it != sparsity.end();
       ++it)
    {
      const auto i = (*it).row();
      const auto j = (*it).column();
      AssertThrow(matrix(i, j) == full(i, j),
                  ExcMessage(std::to_string(matrix(i, j)) +
                             " != " + std::to_string(full(i, j)) + " for i=" +
                             std::to_string(i) + " j=" + std::to_string(j)));
    }

  deallog << "Ok" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  // MPILogInitAll log;
  // test();

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();
      test();
    }
  else
    test();

  return 0;
}
