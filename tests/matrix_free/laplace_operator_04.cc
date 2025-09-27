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



// check the diagonal of LaplaceOperator for the scalar and vector-valued case

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/operators.h>

#include <iostream>

#include "../tests.h"


template <int dim, int fe_degree>
void
test()
{
  using number = double;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_Q<dim>       fe_q(fe_degree);
  FESystem<dim>   fe(fe_q, dim);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  std::shared_ptr<MatrixFree<dim, number>> mf_data(
    new MatrixFree<dim, number>());
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_JxW_values;
    mf_data->reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeOperators::LaplaceOperator<
    dim,
    fe_degree,
    fe_degree + 1,
    1,
    LinearAlgebra::distributed::Vector<number>>
    mf_scalar;
  MatrixFreeOperators::LaplaceOperator<
    dim,
    fe_degree,
    fe_degree + 1,
    dim,
    LinearAlgebra::distributed::Vector<number>>
    mf_vector;
  mf_scalar.initialize(mf_data);
  mf_vector.initialize(mf_data);
  mf_scalar.compute_diagonal();
  mf_vector.compute_diagonal();
  deallog << "Matrix diagonal scalar: " << std::endl;
  mf_scalar.get_matrix_diagonal()->get_vector().print(
    deallog.get_file_stream());
  deallog << "Matrix diagonal vector: " << std::endl;
  mf_vector.get_matrix_diagonal()->get_vector().print(
    deallog.get_file_stream());
}


int
main()
{
  initlog();

  deallog << std::setprecision(4);

  test<2, 2>();
}
