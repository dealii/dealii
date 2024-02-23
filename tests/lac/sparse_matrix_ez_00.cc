// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test SparseMatrixEZ class's behavior on inserting zero entry and the
// elide_zero_values flag.

#include <deal.II/lac/sparse_matrix_ez.h>

#include "../tests.h"

#define NUMBER double

int
main(int argc, char *argv[])
{
  SparseMatrixEZ<NUMBER> m_ez(5, 5, 3);

  std::ofstream logfile("output");
  // Not use because SparseMatrixEZ has no interface with LogStream either with
  // function .print() or operator <<

  // Initialize the matrix
  m_ez.set(0, 0, 1);

  m_ez.set(1, 0, 2);
  m_ez.set(1, 1, 3);

  // The third row is filled to default length
  m_ez.set(2, 0, 7);
  m_ez.set(2, 1, 6.18);
  m_ez.set(2, 2, 3.14);

  // Put a non-zero at end of forth row's default length
  m_ez.set(3, 2, 2.718);

  // Fifth row is left empty

  //------------------------

  logfile << "The initial matrix:" << std::endl;
  m_ez.print(logfile);

  m_ez.set(0, 1, 0.0);
  logfile << "Ignore a zero entry insertion:" << std::endl;
  m_ez.print(logfile);

  m_ez.set(1, 0, 0.0);
  logfile << "Set existing entry (1,0) to zero:" << std::endl;
  m_ez.print(logfile);

  m_ez.set(1, 2, 0.0, /*elide_zero_values=*/false);
  logfile << "Insert new zero to (1,2):" << std::endl;
  m_ez.print(logfile);

  m_ez.set(2, 4, 0.0, /*elide_zero_values=*/false);
  logfile << "Insert new zero to (2,4), out of range of default row length:"
          << std::endl;
  m_ez.print(logfile);

  m_ez.set(2, 1, 0.0, /*elide_zero_values=*/false);
  logfile << "Set existing entry (2,1) to zero with elide_zero_values=false:"
          << std::endl;
  m_ez.print(logfile);

  m_ez.set(3, 0, 1.414, /*elide_zero_values=*/false);
  logfile << "Insert a non-zero entry at (3,0) with elide_zero_values=false:"
          << std::endl;
  m_ez.print(logfile);

  m_ez.set(4, 4, 0.0, /*elide_zero_values=*/false);
  logfile
    << "Insert new zero entry to an empty row with elide_zero_values=false:"
    << std::endl;
  m_ez.print(logfile);

  {
    SparseMatrixEZ<NUMBER> m_ez_cpoier(5, 5, 2);
    m_ez_cpoier.copy_from(m_ez);
    logfile << "Copy the matrix and drop all zeros:" << std::endl;
    m_ez_cpoier.print(logfile);
  }
  {
    SparseMatrixEZ<NUMBER> m_ez_cpoier(5, 5, 2);
    m_ez_cpoier.copy_from(m_ez, /*elide_zero_values=*/false);
    logfile << "Copy the matrix and keep all zeros:" << std::endl;
    m_ez_cpoier.print(logfile);
  }

  logfile.close();
  return (0);
}
