// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// exercise DerivativeForm::norm()


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
      {
        dF[i][j] = (i + 2 * j + 1);
        dF_norm_sqr += (i + 2 * j + 1) * (i + 2 * j + 1);
      }

  DerivativeForm<2, dim, spacedim> ddF;
  double                           ddF_norm_sqr = 0;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        {
          ddF[i][j][k] = (i + 2 * j + 3 * k + 1);
          ddF_norm_sqr += (i + 2 * j + 3 * k + 1) * (i + 2 * j + 3 * k + 1);
        }

  // output the norms of these objects
  deallog << "||dF||: " << dF.norm() << std::endl;
  deallog << "||ddF||: " << ddF.norm() << std::endl;

  Assert(std::fabs(dF.norm() - std::sqrt(dF_norm_sqr)) < 1e-12,
         ExcInternalError());
  Assert(std::fabs(ddF.norm() - std::sqrt(ddF_norm_sqr)) < 1e-12,
         ExcInternalError());
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
