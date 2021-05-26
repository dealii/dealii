// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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
  for (unsigned int i = 1; i < 25; i++)
    {
      for (unsigned int j = 1; j < 25; j++)
        if (FEEvaluation<dim, -1, 0, 1>::fast_evaluation_supported(i, j))
          deallog << 1 << " ";
        else
          deallog << "  ";
      deallog << std::endl;
    }
  for (unsigned int i = 1; i < 25; i++)
    {
      for (unsigned int j = 1; j < 25; j++)
        if (FEFaceEvaluation<dim, -1, 0, 1>::fast_evaluation_supported(i, j))
          deallog << 1 << " ";
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
