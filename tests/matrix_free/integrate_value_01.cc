// ---------------------------------------------------------------------
//
// Copyright (C) 2024 by the deal.II authors
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



// Test FEEvaluation::integrate_value() and
// FEPointEvaluation::integrate_value().

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

template <int dim, typename Number>
class RightHandSideFunction : public Function<dim, Number>
{
public:
  RightHandSideFunction(const unsigned int n_components)
    : Function<dim, Number>(n_components)
  {}

  Number
  value(const Point<dim, Number> &p, const unsigned int component) const
  {
    if (component == 0)
      return p[0];
    else if (component == 1)
      return p[1] * p[1];
    else
      return 0;
  }

  VectorizedArray<Number>
  value(const Point<dim, VectorizedArray<Number>> &p,
        const unsigned int                         component) const
  {
    if (component == 0)
      return p[0];
    else if (component == 1)
      return p[1] * p[1];
    else
      return 0;
  }
};


template <typename Number>
Number
sum(const Number &value)
{
  return value;
}


template <typename Number, int dim>
Number
sum(const Tensor<1, dim, Number> &value)
{
  Number result = {};

  for (unsigned int i = 0; i < dim; ++i)
    result += value[i];

  return result;
}


template <typename Number>
Number &
get(Number &value, const unsigned int)
{
  return value;
}


template <typename Number, int dim>
Number &
get(Tensor<1, dim, Number> &value, const unsigned int component)
{
  return value[component];
}



template <int dim, int n_components>
void
test()
{
  const unsigned int fe_degree                 = 2;
  const unsigned int n_global_mesh_refinements = 2;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_global_mesh_refinements);

  FESystem<dim> fe(FE_Q<dim>(fe_degree), n_components);
  QGauss<dim>   quad(fe_degree + 1);

  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> vector(dof_handler.n_dofs());

  RightHandSideFunction<dim, double> fu(n_components);

  Vector<double> integrals(tria.n_active_cells());

  // create reference solution with VectorTools::integrate_difference()
  {
    Functions::ConstantFunction<dim, double> zero(0.0, n_components);
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      vector,
                                      fu,
                                      integrals,
                                      quad,
                                      VectorTools::NormType::L2_norm);

    for (unsigned int i = 0; i < integrals.size(); ++i)
      integrals[i] = integrals[i] * integrals[i];
  }

  // create solution with FEEvaluation::integrate_value()
  {
    Vector<double> my_integrals(tria.n_active_cells());

    typename MatrixFree<dim>::AdditionalData ad;
    ad.mapping_update_flags =
      update_JxW_values | update_values | update_quadrature_points;
    MatrixFree<dim> matrix_free;

    matrix_free.reinit(
      mapping, dof_handler, AffineConstraints<double>(), quad, ad);

    FEEvaluation<dim, -1, 0, n_components> fe_eval(matrix_free);

    for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
      {
        fe_eval.reinit(cell);

        for (const auto q : fe_eval.quadrature_point_indices())
          {
            const auto point = fe_eval.quadrature_point(q);

            typename FEEvaluation<dim, -1, 0, n_components>::value_type value;

            for (unsigned int c = 0; c < n_components; ++c)
              get(value, c) = fu.value(point, c) * fu.value(point, c);

            fe_eval.submit_value(value, q);
          }

        const auto integrated_cell_value = sum(fe_eval.integrate_value());

        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_cell_batch(cell);
             ++v)
          my_integrals[matrix_free.get_cell_iterator(cell, v)
                         ->active_cell_index()] = integrated_cell_value[v];
      }

    for (unsigned int i = 0; i < my_integrals.size(); ++i)
      Assert(std::abs(my_integrals[i] - integrals[i]) < 1e-12,
             ExcInternalError());
  }

  // create solution with FEEvaluation::integrate_value() (non-vectorized)
  {
    Vector<double> my_integrals(tria.n_active_cells());

    NonMatching::MappingInfo<dim, dim, double> mapping_info(mapping,
                                                            update_JxW_values |
                                                              update_values);

    std::vector<typename Triangulation<dim>::cell_iterator> cells;
    std::vector<Quadrature<dim>>                            quadratures;

    for (const auto &cell : tria.active_cell_iterators())
      {
        cells.push_back(cell);
        quadratures.push_back(quad);
      }

    mapping_info.reinit_cells(cells, quadratures);

    FEPointEvaluation<n_components, dim, dim, double> fe_eval(mapping_info, fe);

    for (unsigned int cell = 0; cell < cells.size(); ++cell)
      {
        fe_eval.reinit(cell);

        for (const auto q : fe_eval.quadrature_point_indices())
          {
            const auto point = fe_eval.real_point(q);

            typename FEPointEvaluation<n_components, dim, dim, double>::
              value_type value;

            for (unsigned int c = 0; c < n_components; ++c)
              get(value, c) = fu.value(point, c) * fu.value(point, c);

            fe_eval.submit_value(value, q);
          }

        my_integrals[cell] = sum(fe_eval.integrate_value());
      }

    for (unsigned int i = 0; i < my_integrals.size(); ++i)
      Assert(std::abs(my_integrals[i] - integrals[i]) < 1e-12,
             ExcInternalError());
  }

  // create solution with FEEvaluation::integrate_value() (vectorized)
  {
    Vector<double> my_integrals(tria.n_active_cells());

    NonMatching::MappingInfo<dim, dim, VectorizedArray<double>> mapping_info(
      mapping, update_JxW_values | update_values);

    std::vector<typename Triangulation<dim>::cell_iterator> cells;
    std::vector<Quadrature<dim>>                            quadratures;

    for (const auto &cell : tria.active_cell_iterators())
      {
        cells.push_back(cell);
        quadratures.push_back(quad);
      }

    mapping_info.reinit_cells(cells, quadratures);

    FEPointEvaluation<n_components, dim, dim, VectorizedArray<double>> fe_eval(
      mapping_info, fe);

    for (unsigned int cell = 0; cell < cells.size(); ++cell)
      {
        fe_eval.reinit(cell);

        for (const auto q : fe_eval.quadrature_point_indices())
          {
            const auto point = fe_eval.real_point(q);

            typename FEPointEvaluation<n_components,
                                       dim,
                                       dim,
                                       VectorizedArray<double>>::value_type
              value;

            for (unsigned int c = 0; c < n_components; ++c)
              get(value, c) = fu.value(point, c) * fu.value(point, c);

            fe_eval.submit_value(value, q);
          }

        my_integrals[cell] = sum(fe_eval.integrate_value());
      }

    for (unsigned int i = 0; i < my_integrals.size(); ++i)
      Assert(std::abs(my_integrals[i] - integrals[i]) < 1e-12,
             ExcInternalError());
  }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2, 1>();
  test<2, 2>();
  test<2, 3>();

  test<3, 1>();
  test<3, 2>();
  test<3, 3>();
  test<3, 4>();
}
