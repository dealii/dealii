// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test copying between matrices (copy_from and fill)

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  SparseMatrixEZ<double> ez(5, 4);
  ez.set(0, 0, 2.);
  ez.set(0, 2, 3.);
  ez.set(0, 3, 4.);
  ez.set(1, 0, 5.);
  ez.set(1, 1, 6.);
  ez.set(1, 3, 7.);
  ez.set(2, 0, 8.);
  ez.set(2, 1, 9.);
  ez.set(2, 2, 10.);
  ez.set(2, 3, 11.);
  ez.set(4, 0, 12.);
  ez.set(4, 2, 13.);
  ez.set(4, 3, 14.);

  deallog << "FullMatrix<float>::copy_from  SparseMatrixEZ<double>"
          << std::endl;
  FullMatrix<float> ff;
  ff.copy_from(ez);
  ff.print_formatted(deallog.get_file_stream(), 0, false, 5, "~");

  deallog << "LAPACKFullMatrix<double>::copy_from  SparseMatrixEZ<double>"
          << std::endl;
  LAPACKFullMatrix<double> lfd;
  lfd.copy_from(ez);
  lfd.print_formatted(deallog.get_file_stream(), 0, false, 5, "~");

  lfd.reinit(2, 3);
  deallog << "LAPACKFullMatrix<double>::fill  SparseMatrixEZ<double>"
          << std::endl;
  lfd.fill(ez, 0, 0, 2, 1);
  lfd.print_formatted(deallog.get_file_stream(), 0, false, 5, "~");
  deallog << "LAPACKFullMatrix<double>::fill  SparseMatrixEZ<double>"
          << std::endl;
  lfd.fill(ez, 1, 1, 4, 2);
  lfd.print_formatted(deallog.get_file_stream(), 0, false, 5, "~");
}
