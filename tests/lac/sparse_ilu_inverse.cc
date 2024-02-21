// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// compute the inverse of a small matrix using the SparseILU with
// infinite fill-in

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);


  for (unsigned int N = 1; N < 5; ++N)
    {
      deallog << "N=" << N << std::endl;

      SparsityPattern structure(N, N, N);
      for (unsigned int i = 0; i < N; ++i)
        for (unsigned int j = 0; j < N; ++j)
          structure.add(i, j);
      structure.compress();
      SparseMatrix<double> A(structure);
      for (unsigned int i = 0; i < N; ++i)
        {
          A.set(i, i, 2);
          if (i >= 1)
            A.set(i, i - 1, -1);
          if (i < N - 1)
            A.set(i, i + 1, -1);
        }

      SparseILU<double> ilu;
      ilu.initialize(A, SparseILU<double>::AdditionalData());

      // now get an explicit
      // representation of the
      // inverse
      FullMatrix<double> inverse(N, N);
      Vector<double>     tmp1(N), tmp2(N);
      for (unsigned int i = 0; i < N; ++i)
        {
          tmp1    = 0;
          tmp1(i) = 1;
          ilu.vmult(tmp2, tmp1);
          for (unsigned int j = 0; j < N; ++j)
            inverse(i, j) = tmp2(j);
        }

      deallog << "Matrix A:" << std::endl;
      A.print_formatted(deallog.get_file_stream(), 3, false);

      deallog << "Matrix A^{-1}:" << std::endl;
      inverse.print_formatted(deallog.get_file_stream(), 3, false);

      deallog << std::endl;
    }
}
