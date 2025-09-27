// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This test template evaluates a simple operator on FE_PolyTensor
// elements using FEFaceEvaluation and compares the result with the output
// of FEFaceValues (which is considered to be the reference) on cartesian
// meshes without hanging nodes. (It will be extended to also handle general
// meshes and hanging nodes in the future.) The test do not include
// multithreading because FEFaceValues is not thread-safe.
// See matrix_vector_rt_face_01.cc for an example that uses this template.

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

enum TestType : unsigned char
{
  values           = 0,
  values_gradients = 1,
  gradients        = 2,
  divergence       = 3
};

std::string
enum_to_string(const TestType enum_type)
{
  std::string string_type;
  switch (enum_type)
    {
      case TestType::values:
        string_type = "Values ";
        break;
      case TestType::gradients:
        string_type = "Gradients ";
        break;
      case TestType::values_gradients:
        string_type = "Values and Gradients ";
        break;
      case TestType::divergence:
        string_type = "Divergence ";
        break;
      default:
        AssertThrow(false, ExcNotImplemented());
        break;
    }
  return string_type;
}

template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in,
                 const TestType                 test_type)
    : data(data_in)
    , test_type(test_type)
  {
    evaluation_flag =
      (test_type == TestType::values) ?
        EvaluationFlags::values :
        ((test_type == TestType::gradients) ?
           EvaluationFlags::gradients :
           ((test_type == TestType::values_gradients) ?
              EvaluationFlags::values | EvaluationFlags::gradients :
              EvaluationFlags::gradients));
  };

  virtual ~MatrixFreeTest(){};

  void
  dummy(const MatrixFree<dim, Number> &,
        Vector<Number> &,
        const Vector<Number> &,
        const std::pair<unsigned int, unsigned int> &) const
  {}

  void
  operator_face(const MatrixFree<dim, Number>               &data,
                Vector<Number>                              &dst,
                const Vector<Number>                        &src,
                const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_eval(data,
                                                                         true);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_eval_n(
      data, false);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval_n.reinit(face);

        fe_eval.gather_evaluate(src, evaluation_flag);
        fe_eval_n.gather_evaluate(src, evaluation_flag);
        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            if (test_type < TestType::gradients)
              {
                fe_eval.submit_value(10. * fe_eval.get_value(q), q);
                fe_eval_n.submit_value(10. * fe_eval_n.get_value(q), q);
              }
            if (test_type == TestType::gradients ||
                test_type == TestType::values_gradients)
              {
                fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
                fe_eval_n.submit_gradient(fe_eval_n.get_gradient(q), q);
              }
            else if (test_type == TestType::divergence)
              {
                fe_eval.submit_divergence(fe_eval.get_divergence(q), q);
                fe_eval_n.submit_divergence(fe_eval_n.get_divergence(q), q);
              }
          }
        fe_eval.integrate_scatter(evaluation_flag, dst);
        fe_eval_n.integrate_scatter(evaluation_flag, dst);
      }
  };

  void
  operator_boundary(
    const MatrixFree<dim, Number>               &data,
    Vector<Number>                              &dst,
    const Vector<Number>                        &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, dim, Number> fe_eval(data,
                                                                         true);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval.gather_evaluate(src, evaluation_flag);

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            if (test_type < TestType::gradients)
              fe_eval.submit_value(10. * fe_eval.get_value(q), q);
            if (test_type == TestType::gradients ||
                test_type == TestType::values_gradients)
              fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
            else if (test_type == TestType::divergence)
              fe_eval.submit_divergence(fe_eval.get_divergence(q), q);
          }

        fe_eval.integrate_scatter(evaluation_flag, dst);
      }
  };


  void
  test_functions(Vector<Number> &dst, const Vector<Number> &src) const
  {
    data.loop(&MatrixFreeTest::dummy,
              &MatrixFreeTest::operator_face,
              &MatrixFreeTest::operator_boundary,
              this,
              dst,
              src);
  };

protected:
  const MatrixFree<dim, Number>   &data;
  EvaluationFlags::EvaluationFlags evaluation_flag;
  const TestType                   test_type;
};



template <int dim, int fe_degree, typename Number>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<double> &constraints,
        const TestType                   test_type)
{
  deallog << "Testing " << enum_to_string(test_type) << std::endl;

  constexpr unsigned int n_q_points = fe_degree + 2;
  constexpr unsigned int m_degree   = fe_degree + 2;

  const MappingQ<dim>     mapping(m_degree);
  MatrixFree<dim, Number> mf_data;
  {
    const QGauss<1>                                  quad(n_q_points);
    typename MatrixFree<dim, Number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, Number>::AdditionalData::none;
    data.mapping_update_flags  = update_piola;
    data.mapping_update_flags_inner_faces = update_contravariant_transformation;
    mf_data.reinit(mapping, dof, constraints, quad, data);
  }

  Vector<Number> solution, initial_condition;

  // create vector with random entries
  mf_data.initialize_dof_vector(solution);
  mf_data.initialize_dof_vector(initial_condition);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      initial_condition[i] = random_value<Number>();
    }
  constraints.distribute(initial_condition);

  MatrixFreeTest<dim, -1, 0, Number> mf(mf_data, test_type);
  mf.test_functions(solution, initial_condition);


  // Evaluation with FEFaceValues
  SparsityPattern        sp;
  SparseMatrix<double>   system_matrix;
  DynamicSparsityPattern dsp(dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern(dof, dsp, constraints, false);
  sp.copy_from(dsp);
  system_matrix.reinit(sp);

  FEFaceValues<dim> fe_val(mapping,
                           dof.get_fe(),
                           QGauss<dim - 1>(n_q_points),
                           update_values | update_gradients |
                             update_JxW_values | update_piola);

  const unsigned int dofs_per_cell = fe_val.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<Tensor<1, dim>> phi_val(dofs_per_cell);
  std::vector<Tensor<2, dim>> phi_grad(dofs_per_cell);
  std::vector<double>         phi_div(dofs_per_cell);

  const FEValuesExtractors::Vector velocities(0);
  // Assemble matrix
  for (const auto &cell : dof.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      {
        fe_val.reinit(cell, face);
        local_matrix = 0;

        for (const auto q : fe_val.quadrature_point_indices())
          {
            const double JxW = fe_val.JxW(q);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                phi_val[i]  = fe_val[velocities].value(i, q);
                phi_grad[i] = fe_val[velocities].gradient(i, q);
                phi_div[i]  = fe_val[velocities].divergence(i, q);
              }
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  if (test_type < TestType::gradients)
                    local_matrix(i, j) += 10. * (phi_val[j] * phi_val[i]) * JxW;
                  if (test_type == TestType::gradients ||
                      test_type == TestType::values_gradients)
                    local_matrix(i, j) +=
                      scalar_product(phi_grad[i], phi_grad[j]) * JxW;
                  else if (test_type == TestType::divergence)
                    local_matrix(i, j) += phi_div[i] * phi_div[j] * JxW;
                }
          }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_dof_indices,
                                               system_matrix);
      }

  Vector<Number> ref(solution.size());

  // Compute reference
  system_matrix.vmult(ref, initial_condition);
  constraints.set_zero(ref);

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
    // test<3, 3>(); //Takes way too long
    deallog.pop();
  }
}
