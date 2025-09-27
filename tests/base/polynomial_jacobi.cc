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

// check jacobi_polynomial_value and jacobi_polynomial_roots

#include <deal.II/base/polynomial.h>

#include "../tests.h"

using namespace Polynomials;


int
main()
{
  initlog();
  deallog.precision(10);

  for (int alpha = 0; alpha < 3; ++alpha)
    for (int beta = 0; beta < 3; ++beta)
      for (unsigned int degree = 0; degree < 40; ++degree)
        {
          deallog << "Jacobi_" << degree << "^(" << alpha << ',' << beta
                  << ") at 0.3: "
                  << jacobi_polynomial_value(degree, alpha, beta, 0.3)
                  << std::endl;

          deallog << "Roots: ";
          std::vector<double> roots =
            jacobi_polynomial_roots<double>(degree, alpha, beta);

          // assert that roots are increasing
          for (unsigned int i = 1; i < roots.size(); ++i)
            Assert(roots[i] > roots[i - 1], ExcInternalError());
          Assert(roots.size() == 0 || roots.front() > -1., ExcInternalError());
          Assert(roots.size() == 0 || roots.back() < 1., ExcInternalError());

          for (unsigned int i = 0; i < roots.size(); ++i)
            deallog << roots[i] << ' ';
          deallog << std::endl << std::endl;
        }
}
