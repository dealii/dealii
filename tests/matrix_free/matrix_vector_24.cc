// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
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
// matrix-vector products for a special case of linear elements where all DoFs
// are subject to constraints

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "matrix_vector_mf.h"



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType = Vector<Number>,
          int n_q_points_1d   = fe_degree + 1>
class MatrixFreeVariant
{
public:
  using vector_t = VectorizedArray<Number>;

  MatrixFreeVariant(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
  {}

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    const std::function<
      void(const MatrixFree<dim, typename VectorType::value_type> &,
           VectorType &,
           const VectorType &,
           const std::pair<unsigned int, unsigned int> &)>
      wrap = helmholtz_operator<dim, fe_degree, VectorType, n_q_points_1d>;
    data.cell_loop(wrap, dst, src, true);
    for (auto i : data.get_constrained_dofs())
      dst(i) += src(i);
  }

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  using number = double;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeVariant<dim,
                    fe_degree,
                    number,
                    LinearAlgebra::distributed::Vector<number>,
                    fe_degree + 1>
                                             mf(mf_data);
  LinearAlgebra::distributed::Vector<number> in(dof.n_dofs()),
    out(dof.n_dofs());
  LinearAlgebra::distributed::Vector<number> out_dist(in);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      const double entry = random_value<double>();
      in(i)              = entry;
    }

  out = in;
  mf.vmult(out, in);

  // assemble sparse matrix with (\nabla v, \nabla u) + (v, 10 * u)
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
                cell_matrix(i, j) += ((fe_values.shape_grad(i, q_point) *
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
  // set matrix entries to constrained rows to 1 for consistency with
  // matrix-free
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      sparse_matrix.set(i, i, 1);

  sparse_matrix.vmult(out_dist, in);
  out -= out_dist;
  const double diff_norm = out.linfty_norm() / out_dist.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<2, 1>();
  test<3, 1>();
}
