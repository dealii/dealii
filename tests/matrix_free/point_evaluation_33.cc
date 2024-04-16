// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check n_active_entries_per_quadrature_batch for vectorized FEPointEvaluation

#include <deal.II/base/function_lib.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const unsigned int degree, const bool is_mapping_q = true)
{
  Triangulation<dim> tria;

  if (dim > 1)
    GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1, 6);
  else
    GridGenerator::subdivided_hyper_cube(tria, 2, 0, 1);

  std::unique_ptr<Mapping<dim>> mapping;
  if (is_mapping_q)
    mapping = std::make_unique<MappingQ<dim>>(degree);
  else
    mapping = std::make_unique<MappingFE<dim>>(FE_Q<dim>(degree));

  deallog << "Mapping of degree " << degree << " with vectorization length "
          << VectorizedArrayType::size() << std::endl;

  QGauss<dim> quadrature(degree + 1);

  std::vector<Point<dim>> unit_points = quadrature.get_points();

  FE_Q<dim>     fe(degree);
  FEValues<dim> fe_values(*mapping,
                          fe,
                          Quadrature<dim>(unit_points),
                          update_values | update_gradients);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  Vector<Number> vector(dof_handler.n_dofs());

  FEPointEvaluation<1, dim, dim, VectorizedArrayType> evaluator(
    *mapping, fe, update_values | update_gradients);

  Tensor<1, dim, Number> exponents;
  exponents[0] = 1.;
  VectorTools::interpolate(*mapping,
                           dof_handler,
                           Functions::Monomial<dim, Number>(exponents),
                           vector);

  std::vector<Number> solution_values(fe.dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      evaluator.reinit(cell, unit_points);

      for (const unsigned int q : evaluator.quadrature_point_indices())
        deallog << evaluator.n_active_entries_per_quadrature_batch(q) << " ";
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<1>(1);
  test<1>(3);
  test<2>(1);
  test<2>(2);
  test<2>(6);
  test<3>(1);
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS > 0
  test<3, double, VectorizedArray<double, 2>>(5);
  test<3, float, VectorizedArray<float, 4>>(5);
#endif
  test<3>(1, false);
}
