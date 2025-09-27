// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests CellwiseInverseMassMatrix::transform_from_q_points_to_basis

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0.2, 0.8);

  const unsigned int degree = fe_degree;
  FE_DGQ<dim>        fe(degree);
  deallog << "Testing " << fe.get_name() << std::endl;
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  Vector<double>                 vec(dof.n_dofs());
  Functions::SquareFunction<dim> function;
  VectorTools::interpolate(dof, function, vec);

  MatrixFree<dim>                          matrix_free;
  AffineConstraints<double>                dummy;
  typename MatrixFree<dim>::AdditionalData additional_data;
  additional_data.mapping_update_flags |= update_quadrature_points;
  matrix_free.reinit(
    MappingQ1<dim>{}, dof, dummy, QGauss<1>(degree + 1), additional_data);

  FEEvaluation<dim, degree> fe_eval(matrix_free);
  FEEvaluation<dim, degree> fe_eval_default(matrix_free);
  MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree> inverse_mass(
    fe_eval);

  fe_eval.reinit(0);
  fe_eval_default.reinit(0);
  fe_eval_default.read_dof_values(vec);
  for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = fe_eval.quadrature_point(q)[d][0];
      fe_eval.begin_values()[q] = function.value(p);
    }
  inverse_mass.transform_from_q_points_to_basis(1,
                                                fe_eval.begin_values(),
                                                fe_eval.begin_dof_values());
  for (unsigned int q = 0; q < fe_eval.dofs_per_cell; ++q)
    deallog << "Value reference / actual: "
            << fe_eval_default.begin_dof_values()[q][0] << " / "
            << fe_eval.begin_dof_values()[q][0] << std::endl;

  deallog << std::endl;
}



int
main()
{
  initlog();

  deallog << std::setprecision(10);

  {
    deallog.push("2d");
    test<2, 2>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 2>();
    test<3, 3>();
    deallog.pop();
  }
}
