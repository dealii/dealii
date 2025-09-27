// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// evaluating the geometry on the fly as in matrix_vector_15, but using
// vector-valued finite elements.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType = Vector<Number>>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const DoFHandler<dim>           &dof_handler,
                 const AffineConstraints<double> &constraints)
    : dof_handler(dof_handler)
    , constraints(constraints)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    VectorType src_cpy = src;
    constraints.distribute(src_cpy);
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> fe_eval(
      dof_handler.get_fe(),
      QGauss<1>(fe_degree + 1),
      update_values | update_gradients | update_JxW_values);
    dst = 0;
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src_cpy);
        fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            fe_eval.submit_value(make_vectorized_array<Number>(10.) *
                                   fe_eval.get_value(q),
                                 q);
            fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
          }
        fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        fe_eval.distribute_local_to_global(dst);
      }
    constraints.condense(dst);
  };

private:
  const DoFHandler<dim>           &dof_handler;
  const AffineConstraints<double> &constraints;
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
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  MatrixFreeTest<dim, fe_degree, number> mf(dof, constraints);
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


  // assemble sparse matrix with vector-valued form of (\nabla v, \nabla u) +
  // (v, 10 * u)
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
                              update_JxW_values);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

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
                if (dof.get_fe().system_to_component_index(i).first ==
                    dof.get_fe().system_to_component_index(j).first)
                  cell_matrix(i, j) +=
                    ((fe_values.shape_grad(i, q_point) *
                        fe_values.shape_grad(j, q_point) +
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
  GridGenerator::hyper_ball(tria);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const unsigned int degree = fe_degree;
  FESystem<dim>      fe(FE_Q<dim>(degree), dim);
  DoFHandler<dim>    dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(dim),
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
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
