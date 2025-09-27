// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// perform some boundary and face operation but fully overwrite the result
// within the post operation

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
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
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class Matrix
{
public:
  Matrix(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
  {}

  void
  do_nothing(const MatrixFree<dim, Number> &,
             LinearAlgebra::distributed::Vector<Number> &,
             const LinearAlgebra::distributed::Vector<Number> &,
             const std::pair<unsigned int, unsigned int> &) const {};

  void
  do_boundary(const MatrixFree<dim, Number>                    &mf,
              LinearAlgebra::distributed::Vector<Number>       &dst,
              const LinearAlgebra::distributed::Vector<Number> &src,
              const std::pair<unsigned int, unsigned int>      &range) const
  {
    using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, 1, Number>;
    FEFaceIntegrator eval(mf, true);

    for (unsigned int face = range.first; face < range.second; ++face)
      {
        eval.reinit(face);
        eval.gather_evaluate(src, EvaluationFlags::values);
        for (unsigned int q = 0; q < eval.n_q_points; ++q)
          {
            eval.submit_value(1.23 * eval.get_value(q), q);
          }
        eval.integrate_scatter(EvaluationFlags::values, dst);
      }
  };

  void
  do_faces(const MatrixFree<dim, Number>                    &mf,
           LinearAlgebra::distributed::Vector<Number>       &dst,
           const LinearAlgebra::distributed::Vector<Number> &src,
           const std::pair<unsigned int, unsigned int>      &range) const
  {
    using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, 1, Number>;
    FEFaceIntegrator eval(mf, true);

    for (unsigned int face = range.first; face < range.second; ++face)
      {
        eval.reinit(face);
        eval.gather_evaluate(src, EvaluationFlags::values);
        for (unsigned int q = 0; q < eval.n_q_points; ++q)
          {
            eval.submit_value(1.23 * eval.get_value(q), q);
          }
        eval.integrate_scatter(EvaluationFlags::values, dst);
      }
  };

  void
  vmult_pre_post_op_only(
    LinearAlgebra::distributed::Vector<Number>       &dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    data.loop(
      &Matrix::do_nothing,
      &Matrix::do_faces,
      &Matrix::do_boundary,
      this,
      dst,
      src,
      // operation before
      [&](const unsigned int start_range, const unsigned int end_range) {
        for (unsigned int i = start_range; i < end_range; ++i)
          {
            dst.local_element(i) = 0;
          }
      },
      // operation after
      [&](const unsigned int start_range, const unsigned int end_range) {
        for (unsigned int i = start_range; i < end_range; ++i)
          {
            dst.local_element(i) = 2 * src.local_element(i);
          }
      });
  }

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree, typename number>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(7 - dim);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.mapping_update_flags_inner_faces = update_JxW_values;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  Matrix<dim, fe_degree, fe_degree + 1, number> mf(mf_data);
  LinearAlgebra::distributed::Vector<number>    vec1, vec2, vec3;
  mf_data.initialize_dof_vector(vec1);
  mf_data.initialize_dof_vector(vec2);

  for (unsigned int i = 0; i < vec1.locally_owned_size(); ++i)
    {
      // Multiply by 0.01 to make float error with roundoff less than the
      // numdiff absolute tolerance
      double entry          = 0.01 * random_value<double>();
      vec1.local_element(i) = entry;
      entry                 = 0.01 * random_value<double>();
      vec2.local_element(i) = entry;
    }

  mf.vmult_pre_post_op_only(vec2, vec1);

  auto ref = vec1;
  ref *= 2.0;
  vec2 -= ref;
  deallog << "Error: " << vec2.linfty_norm() << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  mpi_initlog();
  deallog << std::setprecision(3);

  test<2, 3, double>();
  test<2, 3, float>();
  test<3, 3, double>();
}
