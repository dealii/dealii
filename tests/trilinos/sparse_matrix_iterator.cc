// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  TrilinosWrappers::SparsityPattern pattern(4, 5, 2);
  pattern.add(0, 2);
  pattern.add(0, 0);
  pattern.add(1, 0);
  pattern.add(1, 3);
  pattern.add(2, 4);
  pattern.add(2, 2);
  pattern.add(3, 0);
  pattern.add(3, 4);
  pattern.compress();

  TrilinosWrappers::SparseMatrix matrix(pattern);
  matrix.set(0, 2, 3.5);
  matrix.set(0, 0, 1.);
  matrix.set(1, 0, -2.);
  matrix.set(1, 3, 1.5);
  matrix.set(2, 4, -2.25);
  matrix.set(2, 2, -0.5);
  matrix.set(3, 0, 2.);
  matrix.set(3, 4, 0.);

  // Print the matrix
  for (TrilinosWrappers::SparseMatrix::const_iterator i = matrix.begin();
       i != matrix.end();
       ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  deallog << std::endl;

  // Add 0.5 to each element
  for (TrilinosWrappers::SparseMatrix::iterator i = matrix.begin();
       i != matrix.end();
       ++i)
    i->value() += .5;

  // Print the matrix
  for (TrilinosWrappers::SparseMatrix::const_iterator i = matrix.begin();
       i != matrix.end();
       ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  deallog << std::endl;

  // Subtract 1 from each element in row 2
  for (TrilinosWrappers::SparseMatrix::iterator i = matrix.begin(2);
       i != matrix.end(2);
       ++i)
    i->value() -= 1.;

  //  Double each element in row 1
  for (TrilinosWrappers::SparseMatrix::iterator i = matrix.begin(1);
       i != matrix.end(1);
       ++i)
    i->value() *= 2;

  // Set the first entry to zero
  matrix.begin()->value() = 0;

  // Print the matrix
  for (TrilinosWrappers::SparseMatrix::const_iterator i = matrix.begin();
       i != matrix.end();
       ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  deallog << std::endl;
}
