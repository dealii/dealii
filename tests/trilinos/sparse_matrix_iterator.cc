// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console (0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  TrilinosWrappers::SparsityPattern pattern(4,5,2);
  pattern.add(0,2);
  pattern.add(0,0);
  pattern.add(1,0);
  pattern.add(1,3);
  pattern.add(2,4);
  pattern.add(2,2);
  pattern.add(3,0);
  pattern.add(3,4);
  pattern.compress();

  TrilinosWrappers::SparseMatrix matrix(pattern);
  matrix.set(0,2, 3.5);
  matrix.set(0,0, 1.);
  matrix.set(1,0, -2.);
  matrix.set(1,3, 1.5);
  matrix.set(2,4, -2.25);
  matrix.set(2,2, -0.5);
  matrix.set(3,0, 2.);
  matrix.set(3,4, 0.);

  // Print the matrix
  for (TrilinosWrappers::SparseMatrix::const_iterator i=matrix.begin(); i != matrix.end(); ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  deallog << std::endl;

  // Add 0.5 to each element
  for (TrilinosWrappers::SparseMatrix::iterator i=matrix.begin(); i != matrix.end(); ++i)
    i->value() += .5;

  // Print the matrix
  for (TrilinosWrappers::SparseMatrix::const_iterator i=matrix.begin(); i != matrix.end(); ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  deallog << std::endl;

  // Subtract 1 from each element in row 2
  for (TrilinosWrappers::SparseMatrix::iterator i=matrix.begin(2); i != matrix.end(2); ++i)
    i->value() -= 1.;

  //  Double each element in row 1
  for (TrilinosWrappers::SparseMatrix::iterator i=matrix.begin(1); i != matrix.end(1); ++i)
    i->value() *= 2;

  // Set the first entry to zero
  matrix.begin()->value() = 0;

  // Print the matrix
  for (TrilinosWrappers::SparseMatrix::const_iterator i=matrix.begin(); i != matrix.end(); ++i)
    deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  deallog << std::endl;
}
