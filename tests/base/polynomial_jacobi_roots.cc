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

// check jacobi_polynomial_roots by comparison to eigenvalue problem with
// LAPACK

#include <deal.II/base/polynomial.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include "../tests.h"

using namespace Polynomials;


int
main()
{
  initlog();

  for (int alpha = 0; alpha < 4; ++alpha)
    for (unsigned int degree = 1; degree < 40; degree += 3)
      {
        const int beta = alpha;
        deallog << "Jacobi_" << degree << "^(" << alpha << ',' << beta << ')'
                << std::endl;

        std::vector<double> roots =
          jacobi_polynomial_roots<double>(degree, alpha, beta);
        AssertDimension(roots.size(), degree);
        Vector<double>           roots_reference(degree);
        LAPACKFullMatrix<double> eigenvalue_mat(degree, degree);
        for (unsigned int k = 1; k < degree; ++k)
          {
            eigenvalue_mat(k - 1, k) = std::sqrt(
              4. * k * (k + alpha) * (k + beta) * (k + alpha + beta) /
              ((2. * k - 1 + alpha + beta) * (2. * k + alpha + beta) *
               (2. * k + alpha + beta) * (2. * k + 1 + alpha + beta)));
            eigenvalue_mat(k, k - 1) = eigenvalue_mat(k - 1, k);
            // diagonal entry is zero for alpha=beta
          }
        FullMatrix<double> eigenvectors(roots_reference.size(),
                                        roots_reference.size());
        eigenvalue_mat.compute_eigenvalues_symmetric(
          -1., 1., 1.e-20, roots_reference, eigenvectors);
        deallog << "Roots (implemented vs reference)" << std::endl;
        for (unsigned int i = 0; i < degree; ++i)
          deallog << std::setw(22) << std::setprecision(16) << roots[i]
                  << " (fval = " << std::setw(9) << std::setprecision(3)
                  << jacobi_polynomial_value(degree, alpha, beta, roots[i])
                  << " ) " << std::setw(22) << std::setprecision(16)
                  << 0.5 + 0.5 * roots_reference[i]
                  << " (fval = " << std::setw(9) << std::setprecision(3)
                  << jacobi_polynomial_value(degree,
                                             alpha,
                                             beta,
                                             0.5 + 0.5 * roots_reference[i])
                  << " )" << std::endl;
        deallog << std::endl;
      }
}
