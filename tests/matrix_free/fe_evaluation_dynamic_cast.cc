// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
