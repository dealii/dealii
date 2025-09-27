// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// exercise DerivativeForm::norm() for complex data types


#include <deal.II/base/derivative_form.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  DerivativeForm<1, dim, spacedim, std::complex<double>> dF;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        dF[i][j] = std::complex<double>(i + 2 * j + 1, i + 2 * j + 1);
      }

  // output the determinants of these objects
  deallog << "det(dF): " << dF.determinant().real() << ' '
          << dF.determinant().imag() << std::endl;
}

int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
