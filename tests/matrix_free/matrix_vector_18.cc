// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this tests matrix-free matrix-vector products in 1D (otherwise similar as
// the other matrix-vector tests)


#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType = Vector<Number>>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.cell_loop(
      &MatrixFreeTest<dim, fe_degree, Number, VectorType>::local_operation,
      this,
      dst,
      src);
  };

private:
  const MatrixFree<dim, Number> &data;

  void
  local_operation(const MatrixFree<dim, Number>               &data,
                  VectorType                                  &out,
                  const VectorType                            &in,
                  const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> fe_eval(data);
    const unsigned int n_q_points = fe_eval.n_q_points;

    Tensor<1, dim, VectorizedArray<Number>> ones;
    for (unsigned int d = 0; d < dim; ++d)
      ones[d] = Number(1);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(in);
        fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                         EvaluationFlags::hessians);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            fe_eval.submit_value(Number(10) * fe_eval.get_value(q), q);
            fe_eval.submit_gradient(fe_eval.get_gradient(q) +
                                      make_vectorized_array<Number>(3.2221) *
                                        (fe_eval.get_hessian(q) * ones),
                                    q);
          }
        fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        fe_eval.distribute_local_to_global(out);
      }
  }
};



template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<double> &constraints,
        const unsigned int               parallel_option = 0)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  if (parallel_option > 0)
    deallog << "Parallel option: " << parallel_option << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    if (parallel_option == 1)
      data.tasks_parallel_scheme =
        MatrixFree<dim, number>::AdditionalData::partition_color;
    else if (parallel_option == 2)
      data.tasks_parallel_scheme =
        MatrixFree<dim, number>::AdditionalData::color;
    else
      {
        Assert(parallel_option == 0, ExcInternalError());
        data.tasks_parallel_scheme =
          MatrixFree<dim, number>::AdditionalData::partition_partition;
      }
    data.tasks_block_size = 7;
    data.mapping_update_flags |= update_hessians;

    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, number> mf(mf_data);
  Vector<number>                         in(dof.n_dofs()), out(dof.n_dofs());
  Vector<number>                         in_dist(dof.n_dofs());
  Vector<number>                         out_dist(in_dist);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = random_value<double>();
      in(i)              = entry;
      in_dist(i)         = entry;
    }

  mf.vmult(out_dist, in_dist);


  // assemble sparse matrix with (\nabla v, \nabla u + 3.2221 * \nabla^2 u *
  // ones) + (v, 10 * u)
  SparsityPattern sparsity;
  {
    DynamicSparsityPattern csp(dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern(dof, csp, constraints, true);
    sparsity.copy_from(csp);
  }
  SparseMatrix<double> sparse_matrix(sparsity);
  {
    QGauss<dim> quadrature_formula(fe_degree + 1);

    FEValues<dim> fe_values(dof.get_fe(),
                            quadrature_formula,
                            update_values | update_gradients |
                              update_JxW_values | update_hessians);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    Tensor<1, dim> ones;
    for (unsigned int d = 0; d < dim; ++d)
      ones[d] = 1.;

    typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                   endc = dof.end();
    for (; cell != endc; ++cell)
      {
        cell_matrix = 0;
        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  ((fe_values.shape_grad(i, q_point) *
                      (fe_values.shape_grad(j, q_point) +
                       3.2221 * (fe_values.shape_hessian(j, q_point) * ones)) +
                    10. * fe_values.shape_value(i, q_point) *
                      fe_values.shape_value(j, q_point)) *
                   fe_values.JxW(q_point));
            }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_matrix,
                                               local_dof_indices,
                                               sparse_matrix);
      }
  }

  sparse_matrix.vmult(out, in);
  out -= out_dist;
  const double diff_norm = out.linfty_norm() / out_dist.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  do_test<dim, fe_degree, double>(dof, constraints);
}



int
main()
{
  initlog();

  deallog << std::setprecision(3);

  {
    deallog.push("1d");
    test<1, 1>();
    test<1, 2>();
    test<1, 3>();
    deallog.pop();
    deallog.push("2d");
    test<2, 1>();
    deallog.pop();
  }
}
