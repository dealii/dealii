// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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
  deallog << "det(dF): " << dF.determinant().real() << " "
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
