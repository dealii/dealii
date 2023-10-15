// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// test operator<< for DerivativeForm


#include <deal.II/base/derivative_form.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  DerivativeForm<1, dim, spacedim> dF;
  double                           dF_norm_sqr = 0;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      dF[i][j] = (i + 2 * j + 1);


  DerivativeForm<2, dim, spacedim> ddF;
  double                           ddF_norm_sqr = 0;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        ddF[i][j][k] = (i + 2 * j + 3 * k + 1);


  // output the norms of these objects
  deallog << "dF : " << dF << std::endl;
  deallog << "ddF: " << ddF << std::endl;
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
