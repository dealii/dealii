// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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


// This test template evaluates a simple operator on FE_PolyTensor
// elements using FEEvaluation and compares the result with the output
// of FEValues (which is considered to be the reference) on cartesian
// meshes without hanging nodes. (It will be extended to also handle general
// meshes and hanging nodes in the future.) The test do not include
// multithreading because FEValues is not thread-safe.
// See matrix_vector_rt_01.cc for an example that uses this template.

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/matrix_free.templates.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


// forward declare this function. will be implemented in .cc files
template <int dim, int fe_degree>
void
test();



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  virtual ~MatrixFreeTest(){};

  void
  operator()(const MatrixFree<dim, Number> &              data,
             Vector<Number> &                             dst,
             const Vector<Number> &                       src,
             const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_eval(data);

    // OBS! This will need to be modified once the Piola transform is
    // implemented
    unsigned int n_cells =
      data.get_dof_handler().get_triangulation().n_active_cells();
    Number piola =
      (dim == 2) ? n_cells : Utilities::pow((int)std::cbrt(n_cells), 4);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.gather_evaluate(src,
                                EvaluationFlags::values |
                                  EvaluationFlags::gradients);

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            fe_eval.submit_value(Number(10 * piola) * fe_eval.get_value(q), q);
            fe_eval.submit_gradient(Number(piola) * fe_eval.get_gradient(q), q);
          }

        fe_eval.integrate_scatter(EvaluationFlags::values |
                                    EvaluationFlags::gradients,
                                  dst);
      }
  };


  void
  test_functions(Vector<Number> &dst, const Vector<Number> &src) const
  {
    data.cell_loop(&MatrixFreeTest::operator(), this, dst, src);
  };

protected:
  const MatrixFree<dim, Number> &data;
};


template <int dim, int fe_degree, typename Number>
void
do_test(const DoFHandler<dim> &          dof,
        const AffineConstraints<double> &constraints)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  deallog << "Number of cells: " << dof.get_triangulation().n_active_cells()
          << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl
          << std::endl;


  //   constraints.distribute(solution);
  MatrixFree<dim, Number> mf_data;
  {
    const QGaussLobatto<1>                           quad(fe_degree + 2);
    typename MatrixFree<dim, Number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, Number>::AdditionalData::none;
    data.mapping_update_flags  = update_gradients | update_hessians;
    mf_data.reinit(dof, constraints, quad, data);
  }

  // create vector with random entries
  Vector<Number> solution, initial_condition;

  mf_data.initialize_dof_vector(solution);
  mf_data.initialize_dof_vector(initial_condition);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      initial_condition[i] = random_value<Number>();
    }

  // MatrixFree solution
  MatrixFreeTest<dim, fe_degree, fe_degree + 2, Number> mf(mf_data);
  mf.test_functions(solution, initial_condition);


  SparsityPattern        sp;
  SparseMatrix<double>   system_matrix;
  DynamicSparsityPattern dsp(dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, dsp);
  sp.copy_from(dsp);
  system_matrix.reinit(sp);

  FEValues<dim> fe_val(dof.get_fe(),
                       Quadrature<dim>(mf_data.get_quadrature(0)),
                       update_values | update_gradients | update_JxW_values);


  const unsigned int dofs_per_cell = fe_val.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);


  const FEValuesExtractors::Vector velocities(0);

  // Assemble matrix
  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_val.reinit(cell);
      local_matrix = 0;

      for (const auto q : fe_val.quadrature_point_indices())
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const Tensor<1, dim> phi_i = fe_val[velocities].value(i, q) * 10.;
            const Tensor<2, dim> grad_phi_i = fe_val[velocities].gradient(i, q);

            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const Tensor<1, dim> phi_j = fe_val[velocities].value(j, q);
                const Tensor<2, dim> grad_phi_j =
                  fe_val[velocities].gradient(j, q);

                local_matrix(i, j) +=
                  (phi_j * phi_i + scalar_product(grad_phi_i, grad_phi_j)) *
                  fe_val.JxW(q);
              }
          }
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            local_matrix(i, j));
    }

  Vector<Number> ref(solution.size());

  // Compute reference
  system_matrix.vmult(ref, initial_condition);

  ref -= solution;

  const double diff_norm = ref.linfty_norm() / solution.linfty_norm();
  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2, 2>();
    test<2, 3>();
    deallog.pop();
    deallog.push("3d");
    test<3, 2>();
    test<3, 3>();
    deallog.pop();
  }
}
