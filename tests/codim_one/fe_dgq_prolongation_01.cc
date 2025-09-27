// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this is the root cause for solution_transfer_01: the prolongation
// matrices for FE_DGQ<dim-1,dim> were not computed at all

#include <deal.II/fe/fe_dgq.h>

#include "../tests.h"



int
main()
{
  initlog();

  const unsigned int spacedim = 2;
  const unsigned int dim      = spacedim - 1;

  for (unsigned int degree = 0; degree < 3; ++degree)
    {
      deallog << "Degree=" << degree << std::endl;
      FE_DGQ<dim, spacedim> fe(degree);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
          deallog << fe.get_prolongation_matrix(0)(i, j) << std::endl;
    }

  return 0;
}
