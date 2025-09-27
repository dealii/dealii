// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  AffineConstraints<double> constraints;
  unsigned int              IDs[]  = {1, 2, 3, 5, 8, 13, 21};
  double                    vals[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  for (unsigned int i = 0; i < sizeof(IDs) / sizeof(IDs[0]); ++i)
    {
      constraints.constrain_dof_to_zero(IDs[i]);
      constraints.set_inhomogeneity(IDs[i], vals[i]);
    }

  constraints.print(deallog.get_file_stream());
  deallog << std::endl;

  AffineConstraints<double> cm(std::move(constraints));
  cm.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}
