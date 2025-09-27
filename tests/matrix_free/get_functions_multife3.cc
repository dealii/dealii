// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
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
// operations in getting the function values, the function gradients, and the
// function Laplacians on a hypecube mesh with adaptive refinement with
// components in two different DoFHandlers, one of which is in hp-mode.

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <utility>
#include <vector>

#include "../tests.h"


template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTest
{
public:
  using VectorType = std::vector<Vector<Number>>;

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in,
                 const hp::FECollection<dim>   &fe_collection)
    : data(data_in)
    , fe_val0(data.get_dof_handler(0).get_fe(),
              Quadrature<dim>(data.get_quadrature(0)),
              update_values | update_gradients | update_hessians)
    , hp_fe_val1(fe_collection,
                 hp::QCollection<dim>(data.get_quadrature(1)),
                 update_values | update_gradients | update_hessians){};

  void
  operator()(const MatrixFree<dim, Number> &data,
             VectorType &,
             const VectorType                            &src,
             const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    const auto active_fe_index = data.get_cell_range_category(cell_range);

    FEEvaluation<dim, -1, 0, 1, Number> fe_eval0(data, 0, 0);
    FEEvaluation<dim, -1, 0, 1, Number> fe_eval1(
      data, 1, 1, 0, active_fe_index);
    std::vector<double>         reference_values0(fe_eval0.n_q_points);
    std::vector<Tensor<1, dim>> reference_grads0(fe_eval0.n_q_points);
    std::vector<Tensor<2, dim>> reference_hess0(fe_eval0.n_q_points);
    std::vector<double>         reference_values1(fe_eval1.n_q_points);
    std::vector<Tensor<1, dim>> reference_grads1(fe_eval1.n_q_points);
    std::vector<Tensor<2, dim>> reference_hess1(fe_eval1.n_q_points);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval0.reinit(cell);
        fe_eval0.read_dof_values(src[0]);
        fe_eval0.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                          EvaluationFlags::hessians);

        fe_eval1.reinit(cell);
        fe_eval1.read_dof_values(src[1]);
        fe_eval1.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                          EvaluationFlags::hessians);

        // compare values with the ones the FEValues
        // gives us. Those are seen as reference
        for (unsigned int j = 0; j < data.n_active_entries_per_cell_batch(cell);
             ++j)
          {
            // FE 0
            fe_val0.reinit(data.get_cell_iterator(cell, j, 0));
            fe_val0.get_function_values(src[0], reference_values0);
            fe_val0.get_function_gradients(src[0], reference_grads0);
            fe_val0.get_function_hessians(src[0], reference_hess0);

            for (int q = 0; q < (int)fe_eval0.n_q_points; ++q)
              {
                errors[0] +=
                  std::fabs(fe_eval0.get_value(q)[j] - reference_values0[q]);
                for (unsigned int d = 0; d < dim; ++d)
                  errors[1] += std::fabs(fe_eval0.get_gradient(q)[d][j] -
                                         reference_grads0[q][d]);
                errors[2] += std::fabs(fe_eval0.get_laplacian(q)[j] -
                                       trace(reference_hess0[q]));
                total[0] += std::fabs(reference_values0[q]);
                for (unsigned int d = 0; d < dim; ++d)
                  total[1] += std::fabs(reference_grads0[q][d]);
                total[2] += std::fabs(fe_eval0.get_laplacian(q)[j]);
              }

            // FE 1
            hp_fe_val1.reinit(data.get_cell_iterator(cell, j, 1));
            const auto &fe_val1 = hp_fe_val1.get_present_fe_values();
            fe_val1.get_function_values(src[1], reference_values1);
            fe_val1.get_function_gradients(src[1], reference_grads1);
            fe_val1.get_function_hessians(src[1], reference_hess1);

            for (int q = 0; q < (int)fe_eval1.n_q_points; ++q)
              {
                errors[3] +=
                  std::fabs(fe_eval1.get_value(q)[j] - reference_values1[q]);
                for (unsigned int d = 0; d < dim; ++d)
                  errors[4] += std::fabs(fe_eval1.get_gradient(q)[d][j] -
                                         reference_grads1[q][d]);
                errors[5] += std::fabs(fe_eval1.get_laplacian(q)[j] -
                                       trace(reference_hess1[q]));
                total[3] += std::fabs(reference_values1[q]);
                for (unsigned int d = 0; d < dim; ++d)
                  total[4] += std::fabs(reference_grads1[q][d]);
                total[5] += std::fabs(fe_eval1.get_laplacian(q)[j]);
              }
          }
      }
  }

  void
  test_functions(const VectorType &src) const
  {
    for (unsigned int i = 0; i < 3 * 2; ++i)
      {
        errors[i] = 0;
        total[i]  = 0;
      }
    VectorType dst_dummy;
    data.cell_loop(
      &MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number>::operator(),
      this,
      dst_dummy,
      src);

    // for doubles, use a stricter condition then
    // for floats for the relative error size
    for (unsigned int i = 0; i < 2; ++i)
      {
        if (std::is_same_v<Number, double> == true)
          {
            deallog << "Error function values FE " << i << ": "
                    << errors[i * 3 + 0] / total[i * 3 + 0] << std::endl;
            deallog << "Error function gradients FE " << i << ": "
                    << errors[i * 3 + 1] / total[i * 3 + 1] << std::endl;

            // need to set quite a loose tolerance because
            // FEValues approximates Hessians with finite
            // differences, which are not so
            // accurate. moreover, Hessians are quite
            // large since we chose random numbers. for
            // some elements, it might also be zero
            // (linear elements on quadrilaterals), so
            // need to check for division by 0, too.
            const double output2 =
              total[i * 3 + 2] == 0 ? 0. : errors[i * 3 + 2] / total[i * 3 + 2];
            deallog << "Error function Laplacians FE " << i << ": " << output2
                    << std::endl;
          }
        else if (std::is_same_v<Number, float> == true)
          {
            deallog << "Error function values FE " << i << ": "
                    << errors[i * 3 + 0] / total[i * 3 + 0] << std::endl;
            deallog << "Error function gradients FE " << i << ": "
                    << errors[i * 3 + 1] / total[i * 3 + 1] << std::endl;
            const double output2 =
              total[i * 3 + 2] == 0 ? 0. : errors[i * 3 + 2] / total[i * 3 + 2];
            deallog << "Error function Laplacians FE " << i << ": " << output2
                    << std::endl;
          }
      }
  }

private:
  const MatrixFree<dim, Number> &data;
  mutable FEValues<dim>          fe_val0;
  mutable hp::FEValues<dim>      hp_fe_val1;
  mutable double                 errors[6], total[6];
};



template <int dim, int fe_degree>
void
test()
{
  using number = double;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.begin_active(tria.n_levels() - 2)->set_refine_flag();
  tria.begin_active(tria.n_levels() - 3)->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim> fe0(fe_degree);
  FE_Q<dim> fe1(fe_degree + 1);

  DoFHandler<dim> dof0(tria);
  dof0.distribute_dofs(fe0);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe0);
  fe_collection.push_back(fe1);

  DoFHandler<dim> dof1(tria);
  for (const auto &cell : dof1.active_cell_iterators())
    cell->set_active_fe_index(Testing::rand() % fe_collection.size());
  dof1.distribute_dofs(fe_collection);

  std::vector<const DoFHandler<dim> *> dof(2);
  dof[0] = &dof0;
  dof[1] = &dof1;

  deallog << "Testing " << fe0.get_name() << " and " << fe1.get_name()
          << std::endl;
  // std::cout << "Number of cells: " << tria.n_active_cells() << std::endl;

  std::vector<Vector<double>> src(dof.size());
  for (unsigned int no = 0; no < dof.size(); ++no)
    src[no].reinit(dof[no]->n_dofs());

  std::vector<const AffineConstraints<double> *> constraints(2);
  AffineConstraints<double>                      constraint0;
  DoFTools::make_hanging_node_constraints(*dof[0], constraint0);
  constraint0.close();
  constraints[0] = &constraint0;
  AffineConstraints<double> constraint1;
  DoFTools::make_hanging_node_constraints(*dof[1], constraint1);
  constraint1.close();
  constraints[1] = &constraint1;

  // create vector with random entries
  for (unsigned int no = 0; no < 2; ++no)
    for (unsigned int i = 0; i < dof[no]->n_dofs(); ++i)
      {
        if (constraints[no]->is_constrained(i))
          continue;
        const double entry = random_value<double>();
        src[no](i)         = entry;
      }


  constraints[0]->distribute(src[0]);
  constraints[1]->distribute(src[1]);
  MatrixFree<dim, number> mf_data;
  {
    std::vector<Quadrature<dim>> quad;
    for (unsigned int no = 0; no < 2; ++no)
      quad.push_back(QGauss<dim>(fe_degree + 1 + no));
    mf_data.reinit(MappingQ1<dim>{},
                   dof,
                   constraints,
                   quad,
                   typename MatrixFree<dim, number>::AdditionalData(
                     MatrixFree<dim, number>::AdditionalData::none));
  }

  MatrixFreeTest<dim, fe_degree, fe_degree + 1, number> mf(mf_data,
                                                           fe_collection);
  mf.test_functions(src);
  deallog << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(7);

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
