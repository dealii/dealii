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


// Tests generalized eigenvalues of FullMatrix

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <float.h>

#include <iostream>
#include <vector>

#include "../tests.h"

const double left[4][4] = {{4., -1., -1., -1.},
                           {-1., 4., -1., -1.},
                           {-1., -1., 4., -1.},
                           {-1., -1., -1., 4.}};

const double right[4][4] = {{4., -1., -1., -1.},
                            {-1., 5., -1., -1.},
                            {-1., -1., 6., -1.},
                            {-1., -1., -1., 7.}};



int
main()
{
  initlog();
  deallog.get_file_stream().precision(1);

  FullMatrix<double>       A(4, 4, &left[0][0]), B(4, 4, &right[0][0]);
  LAPACKFullMatrix<double> LA(4, 4), LB(4, 4);
  for (unsigned int itype = 1; itype <= 3; ++itype)
    {
      deallog << std::endl
              << "generalized eigenvalue problem of type " << itype
              << std::endl;
      LA = A;
      LB = B;
      std::vector<Vector<double>> eigenvectors(0);
      Vector<double>              eigenvalues(0);

      LA.compute_generalized_eigenvalues_symmetric(
        LB, 0.5, 3.0, 2.0 * DBL_MIN, eigenvalues, eigenvectors, itype);

      for (unsigned int i = 0; i < eigenvectors.size(); ++i)
        {
          deallog << "generalized eigenvalue " << std::scientific
                  << eigenvalues(i) << std::endl
                  << "generalized eigenvector ";
          for (unsigned int j = 0; j < A.m(); ++j)
            {
              deallog << std::scientific
                      << eigenvectors[i](j) / eigenvectors[i](0) << '\t';
            }
          deallog << std::endl;
        }
    }
}
