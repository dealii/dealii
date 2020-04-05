// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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



// this tests the correctness of matrix free matrix-vector products for two
// vectors on the same DoFHandler. Similar to matrix_vector_12.cc but using
// BlockVector instead of std::vector<Vector>.

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"


template <int dim, int fe_degree, typename Number>
void
helmholtz_operator(const MatrixFree<dim, Number> &                        data,
                   LinearAlgebra::distributed::BlockVector<Number> &      dst,
                   const LinearAlgebra::distributed::BlockVector<Number> &src,
                   const std::pair<unsigned int, unsigned int> &cell_range)
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi0(data);
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi1(data);
  const unsigned int n_q_points = phi0.n_q_points;

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi0.reinit(cell);
      phi1.reinit(cell);

      phi0.read_dof_values(src, 0);
      phi1.read_dof_values(src, 1);
      phi0.evaluate(true, true, false);
      phi1.evaluate(true, true, false);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          phi0.submit_value(make_vectorized_array(Number(10)) *
                              phi0.get_value(q),
                            q);
          phi0.submit_gradient(phi0.get_gradient(q), q);
          phi1.submit_value(make_vectorized_array(Number(10)) *
                              phi1.get_value(q),
                            q);
          phi1.submit_gradient(phi1.get_gradient(q), q);
        }
      phi0.integrate(true, true);
      phi0.distribute_local_to_global(dst, 0);
      phi1.integrate(true, true);
      phi1.distribute_local_to_global(dst, 1);
    }
}



template <int dim, int fe_degree, typename Number>
class MatrixFreeTest
{
public:
  typedef VectorizedArray<Number> vector_t;
  static const std::size_t        n_vectors = VectorizedArray<Number>::size();

  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  void
  vmult(LinearAlgebra::distributed::BlockVector<Number> &      dst,
        const LinearAlgebra::distributed::BlockVector<Number> &src) const
  {
    for (unsigned int i = 0; i < dst.size(); ++i)
      dst[i] = 0;
    const std::function<
      void(const MatrixFree<dim, Number> &,
           LinearAlgebra::distributed::BlockVector<Number> &,
           const LinearAlgebra::distributed::BlockVector<Number> &,
           const std::pair<unsigned int, unsigned int> &)>
      wrap = helmholtz_operator<dim, fe_degree, Number>;
    data.cell_loop(wrap, dst, src);
  };

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree>
void
test()
{
  typedef double number;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  cell                                                   = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->is_locally_owned())
      if (cell->center().norm() < 0.2)
        cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (tria.begin(tria.n_levels() - 1)->is_locally_owned())
    tria.begin(tria.n_levels() - 1)->set_refine_flag();
  if (tria.last()->is_locally_owned())
    tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (cell->is_locally_owned())
          if (counter % (7 - i) == 0)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

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

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  // std::cout << "Number of cells: " << tria.n_global_active_cells() <<
  // std::endl; std::cout << "Number of degrees of freedom: " << dof.n_dofs() <<
  // std::endl; std::cout << "Number of constraints: " <<
  // constraints.n_constraints() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::partition_color;
    data.tasks_block_size = 7;
    mf_data.reinit(dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, number>          mf(mf_data);
  LinearAlgebra::distributed::Vector<number>      ref;
  LinearAlgebra::distributed::BlockVector<number> in(2), out(2);
  for (unsigned int i = 0; i < 2; ++i)
    {
      mf_data.initialize_dof_vector(in.block(i));
      mf_data.initialize_dof_vector(out.block(i));
    }
  in.collect_sizes();
  out.collect_sizes();
  mf_data.initialize_dof_vector(ref);

  for (unsigned int i = 0; i < in.block(0).local_size(); ++i)
    {
      if (constraints.is_constrained(
            dof.locally_owned_dofs().index_within_set(i)))
        continue;
      in.block(0).local_element(i) = random_value<double>();
      in.block(1).local_element(i) = random_value<double>();
    }

  mf.vmult(out, in);


  // assemble sparse matrix with (\nabla v, \nabla u) + (v, 10 * u) for
  // reference
  SparsityPattern sparsity;
  {
    DynamicSparsityPattern dsp(dof.locally_owned_dofs());
    DoFTools::make_sparsity_pattern(dof, dsp, constraints, true);
    sparsity.copy_from(dsp);
  }
  SparseMatrix<number> sparse_matrix(sparsity);
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
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          fe_values.reinit(cell);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
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

  deallog << "Norm of difference (component 1/2): ";
  for (unsigned int i = 0; i < 2; ++i)
    {
      sparse_matrix.vmult(ref, in.block(i));
      out.block(i) -= ref;
      const double diff_norm = out.block(i).linfty_norm();
      deallog << diff_norm << " ";
    }
  deallog << std::endl << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();
  deallog << std::setprecision(4);
  deallog.depth_console(0);

  deallog.push("2d");
  test<2, 1>();
  deallog.pop();

  deallog.push("3d");
  test<3, 1>();
  deallog.pop();
}
