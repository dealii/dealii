// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that FEEvaluationData is virtual.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim    = 2;
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(FE_Q<dim>(1));

  MatrixFree<dim, Number, VectorizedArrayType> dummy;

  dummy.reinit(MappingQ1<dim>{},
               dof,
               AffineConstraints<double>(),
               QGauss<1>(2));

  std::unique_ptr<FEEvaluationData<dim, VectorizedArrayType, false>> phi_base =
    std::make_unique<FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>>(
      dummy);

  if (dynamic_cast<FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> *>(
        phi_base.get()))
    deallog << "OK!" << std::endl;
  else
    deallog << "Fail!" << std::endl;

  if (dynamic_cast<FEEvaluation<dim, -1, 0, 2, Number, VectorizedArrayType> *>(
        phi_base.get()))
    deallog << "Fail!" << std::endl;
  else
    deallog << "OK!" << std::endl;
}
