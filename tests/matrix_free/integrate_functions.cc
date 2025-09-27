// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this function tests the correctness of the implementation of matrix free
// operations in integrating functions and gradients on a hypeball mesh with
// adaptive refinement.

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim, int fe_degree, typename Number>
class MatrixFreeTest
{
public:
  using VectorType = std::vector<Vector<Number> *>;

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
    , fe_val(data.get_dof_handler().get_fe(),
             Quadrature<dim>(data.get_quadrature(0)),
             update_values | update_gradients | update_JxW_values){};

  void
  operator()(const MatrixFree<dim, Number>               &data,
             VectorType                                  &dst,
             const VectorType                            &src,
             const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  test_functions(Vector<Number> &dst, Vector<Number> &dst_deal) const
  {
    dst      = 0;
    dst_deal = 0;
    VectorType dst_data(2);
    dst_data[0] = &dst;
    dst_data[1] = &dst_deal;
    VectorType src_dummy;
    data.cell_loop(&MatrixFreeTest<dim, fe_degree, Number>::operator(),
                   this,
                   dst_data,
                   src_dummy);
  };

private:
  const MatrixFree<dim, Number> &data;
  mutable FEValues<dim>          fe_val;
};



template <int dim, int fe_degree, typename Number>
void
MatrixFreeTest<dim, fe_degree, Number>::operator()(
  const MatrixFree<dim, Number> &data,
  std::vector<Vector<Number> *> &dst,
  const std::vector<Vector<Number> *> &,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> fe_eval(data);
  const unsigned int                     n_q_points    = fe_eval.n_q_points;
  const unsigned int                     dofs_per_cell = fe_eval.dofs_per_cell;
  AlignedVector<VectorizedArray<Number>> values(n_q_points);
  AlignedVector<VectorizedArray<Number>> gradients(dim * n_q_points);
  std::vector<types::global_dof_index>   dof_indices(dofs_per_cell);
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      fe_eval.reinit(cell);
      // compare values with the ones the FEValues
      // gives us. Those are seen as reference
      for (unsigned int j = 0; j < data.n_active_entries_per_cell_batch(cell);
           ++j)
        {
          // generate random numbers at quadrature
          // points and test them with basis functions
          // and their gradients
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              values[q][j] = random_value<double>();
              for (unsigned int d = 0; d < dim; ++d)
                gradients[q * dim + d][j] = -1. + 2. * (random_value<double>());
            }
          fe_val.reinit(data.get_cell_iterator(cell, j));
          data.get_cell_iterator(cell, j)->get_dof_indices(dof_indices);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              double sum = 0.;
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  sum +=
                    values[q][j] * fe_val.shape_value(i, q) * fe_val.JxW(q);
                  for (unsigned int d = 0; d < dim; ++d)
                    sum += (gradients[q * dim + d][j] *
                            fe_val.shape_grad(i, q)[d] * fe_val.JxW(q));
                }
              (*dst[1])(dof_indices[i]) += sum;
            }
        }
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          fe_eval.submit_value(values[q], q);
          Tensor<1, dim, VectorizedArray<Number>> submit;
          for (unsigned int d = 0; d < dim; ++d)
            submit[d] = gradients[q * dim + d];
          fe_eval.submit_gradient(submit, q);
        }
      fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
      fe_eval.distribute_local_to_global(*dst[0]);
    }
}



template <int dim, int fe_degree>
void
test()
{
  using number = double;
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center().norm() < 0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (unsigned int i = 0; i < 7 - 2 * dim; ++i)
    {
      unsigned int counter = 0;
      for (const auto &cell : tria.active_cell_iterators())
        {
          if (counter % (7 - i) == 0)
            cell->set_refine_flag();
          ++counter;
        }
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  deallog << "Testing " << fe.get_name() << std::endl;
  // std::cout << "Number of cells: " << tria.n_active_cells() << std::endl;
  // std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1> quad(fe_degree + 1);
    mf_data.reinit(MappingQ1<dim>{},
                   dof,
                   constraints,
                   quad,
                   typename MatrixFree<dim, number>::AdditionalData(
                     MatrixFree<dim, number>::AdditionalData::none));
  }

  MatrixFreeTest<dim, fe_degree, number> mf(mf_data);
  Vector<number>                         solution(dof.n_dofs());
  Vector<number>                         solution_dist(dof.n_dofs());

  mf.test_functions(solution_dist, solution);

  constraints.condense(solution);

  Vector<number> compare(solution_dist);
  compare -= solution;
  const double diff_norm = compare.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    test<3, 3>();
    deallog.pop();
  }
}
