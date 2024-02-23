// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test FEEvaluation::fast_evaluation_supported()

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"



template <int dim>
void
test()
{
  for (unsigned int i = 1; i < 25; ++i)
    {
      for (unsigned int j = 1; j < 25; ++j)
        if (FEEvaluation<dim, -1, 0, 1>::fast_evaluation_supported(i, j))
          deallog << 1 << ' ';
        else
          deallog << "  ";
      deallog << std::endl;
    }
  for (unsigned int i = 1; i < 25; ++i)
    {
      for (unsigned int j = 1; j < 25; ++j)
        if (FEFaceEvaluation<dim, -1, 0, 1>::fast_evaluation_supported(i, j))
          deallog << 1 << ' ';
        else
          deallog << "  ";
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  test<1>();
}
