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



// test the correctness of the matrix-free cell loop with additional
// operation_before_loop and operation_after_loop lambdas. similar to
// pre_and_post_loops_01 but using Dirichlet boundary conditions as
// constraints

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
  vmult(LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    const std::function<void(const MatrixFree<dim, Number> &,
                             LinearAlgebra::distributed::Vector<Number> &,
                             const LinearAlgebra::distributed::Vector<Number> &,
                             const std::pair<unsigned int, unsigned int> &)>
      wrap = helmholtz_operator<dim,
                                fe_degree,
                                LinearAlgebra::distributed::Vector<Number>,
                                n_q_points_1d>;
    dst    = 0;
    data.cell_loop(wrap, dst, src);
    for (auto i : data.get_constrained_dofs())
      dst.local_element(i) = src.local_element(i);
  }

  void
  vmult_with_update_basic(
    LinearAlgebra::distributed::Vector<Number> &dst,
    LinearAlgebra::distributed::Vector<Number> &src,
    LinearAlgebra::distributed::Vector<Number> &other_vector) const
  {
    src.add(1.5, other_vector);
    other_vector.add(0.7, dst);
    dst = 0;
    vmult(dst, src);
    dst.sadd(0.63, -1.3, other_vector);
  }

  void
  vmult_with_update_merged(
    LinearAlgebra::distributed::Vector<Number> &dst,
    LinearAlgebra::distributed::Vector<Number> &src,
    LinearAlgebra::distributed::Vector<Number> &other_vector) const
  {
    const std::function<void(const MatrixFree<dim, Number> &,
                             LinearAlgebra::distributed::Vector<Number> &,
                             const LinearAlgebra::distributed::Vector<Number> &,
                             const std::pair<unsigned int, unsigned int> &)>
      wrap = helmholtz_operator<dim,
                                fe_degree,
                                LinearAlgebra::distributed::Vector<Number>,
                                n_q_points_1d>;
    data.cell_loop(
      wrap,
      dst,
      src,
      // operation before cell operation
      [&](const unsigned int start_range, const unsigned int end_range) {
        for (unsigned int i = start_range; i < end_range; ++i)
          {
            src.local_element(i) += 1.5 * other_vector.local_element(i);
            other_vector.local_element(i) += 0.7 * dst.local_element(i);
            dst.local_element(i) = 0;
          }
      },
      // operation after cell operation
      [&](const unsigned int start_range, const unsigned int end_range) {
        for (unsigned int i = start_range; i < end_range; ++i)
          {
            dst.local_element(i) =
              0.63 * dst.local_element(i) - 1.3 * other_vector.local_element(i);
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
  const IndexSet relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof);
  constraints.reinit(dof.locally_owned_dofs(), relevant_dofs);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  Matrix<dim, fe_degree, fe_degree + 1, number> mf(mf_data);
  LinearAlgebra::distributed::Vector<number>    vec1, vec2, vec3;
  mf_data.initialize_dof_vector(vec1);
  mf_data.initialize_dof_vector(vec2);
  mf_data.initialize_dof_vector(vec3);

  for (unsigned int i = 0; i < vec1.locally_owned_size(); ++i)
    {
      if (constraints.is_constrained(
            vec1.get_partitioner()->local_to_global(i)))
        continue;

      // Multiply by 0.01 to make float error with roundoff less than the
      // numdiff absolute tolerance
      double entry          = 0.01 * random_value<double>();
      vec1.local_element(i) = entry;
      entry                 = 0.01 * random_value<double>();
      vec2.local_element(i) = entry;
      entry                 = 0.01 * random_value<double>();
      vec3.local_element(i) = entry;
    }

  LinearAlgebra::distributed::Vector<number> ref1 = vec1;
  LinearAlgebra::distributed::Vector<number> ref2 = vec2;
  LinearAlgebra::distributed::Vector<number> ref3 = vec3;

  mf.vmult_with_update_basic(ref3, ref2, ref1);
  mf.vmult_with_update_merged(vec3, vec2, vec1);

  vec3 -= ref3;
  deallog << "Error in 1x merged operation: " << vec3.linfty_norm()
          << std::endl;

  ref3 = 0;

  mf.vmult_with_update_basic(ref1, ref2, ref3);
  mf.vmult_with_update_merged(vec1, vec2, vec3);

  mf.vmult_with_update_basic(ref1, ref2, ref3);
  mf.vmult_with_update_merged(vec1, vec2, vec3);

  vec3 -= ref3;
  deallog << "Error in 2x merged operation: " << vec3.linfty_norm()
          << std::endl;
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
