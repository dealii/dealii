// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Tests eigenvalues of FullMatrix

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <float.h>

#include <iostream>
#include <vector>

#include "../tests.h"

const double left[4][4] = {{1.75, -0.433012701892219, 0.0, 0.0},
                           {-0.433012701892219, 1.25, 0.0, 0.0},
                           {0.0, 0.0, 3.5, -0.5},
                           {0.0, 0.0, -0.5, 3.5}};



int
main()
{
  initlog();
  deallog.get_file_stream().precision(1);

  FullMatrix<double>       A(4, 4, &left[0][0]);
  LAPACKFullMatrix<double> LA(4, 4);
  LA = A;
  FullMatrix<double> eigenvectors;
  Vector<double>     eigenvalues(0);

  LA.compute_eigenvalues_symmetric(
    0.5, 2.5, 2.0 * DBL_MIN, eigenvalues, eigenvectors);

  for (unsigned int i = 0; i < eigenvalues.size(); ++i)
    {
      deallog << "eigenvalue " << std::scientific << eigenvalues(i) << std::endl
              << "eigenvector ";
      for (unsigned int j = 0; j < A.m(); ++j)
        {
          deallog << std::scientific << eigenvectors(j, i) / eigenvectors(0, i)
                  << '\t';
        }
      deallog << std::endl;
    }
}
