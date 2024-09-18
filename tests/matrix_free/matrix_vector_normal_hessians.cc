// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test if applying a matrix vector product with a matrix containing hessians
// on faces produces the same result with FEFaceEvaluation and FEFaceValues.
// This is checked for different combinations of EvaluationFlags, FE types and
// polynomial degrees.

#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test_hessians(const unsigned int                             degree,
              const dealii::FE_Poly<dim>                    &fe,
              const dealii::EvaluationFlags::EvaluationFlags evaluation_flags,
              const bool                                     check_system)
{
  using VectorizedArrayType = VectorizedArray<double>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(4 - dim);

  unsigned int counter = 0;
  for (auto cell = tria.begin_active(); cell != tria.end(); ++cell, ++counter)
    if (cell->is_locally_owned() && counter % 3 == 0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dof_handler(tria);
  if (check_system)
    {
      FESystem<dim> fe_system(fe, dim + 1);
      dof_handler.distribute_dofs(fe_system);
    }
  else
    dof_handler.distribute_dofs(fe);

  MappingQGeneric<dim> mapping(1);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(
    mapping,
    dof_handler,
    0,
    Functions::ZeroFunction<dim>(
      dof_handler.get_fe_collection().n_components()),
    constraints);
  constraints.close();

  // FEFaceEvaluation
  typename MatrixFree<dim, double, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients | update_hessians;

  MatrixFree<dim, double, VectorizedArrayType> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, QGauss<1>(degree + 1), additional_data);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Working with " << fe.get_name() << " and "
            << dof_handler.n_dofs() << " dofs" << std::endl;

  LinearAlgebra::distributed::Vector<double> src, dst, dst2;
  matrix_free.initialize_dof_vector(src);
  for (auto &v : src)
    v = random_value<double>();

  matrix_free.initialize_dof_vector(dst);

  src.update_ghost_values();

  matrix_free.template loop<LinearAlgebra::distributed::Vector<double>,
                            LinearAlgebra::distributed::Vector<double>>(
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval_m(
        matrix_free, true);
      FEFaceEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval_p(
        matrix_free, false);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          // FEFaceEvaluation
          fe_eval_m.reinit(face);
          fe_eval_m.gather_evaluate(src, evaluation_flags);
          for (unsigned int q = 0; q < fe_eval_m.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval_m.submit_normal_hessian(fe_eval_m.get_normal_hessian(q),
                                                q);
            }
          fe_eval_m.integrate(evaluation_flags);
          fe_eval_m.distribute_local_to_global(dst);

          fe_eval_p.reinit(face);
          fe_eval_p.gather_evaluate(src, evaluation_flags);
          for (unsigned int q = 0; q < fe_eval_p.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval_p.submit_normal_hessian(fe_eval_p.get_normal_hessian(q),
                                                q);
            }
          fe_eval_p.integrate(evaluation_flags);
          fe_eval_p.distribute_local_to_global(dst);
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    dst,
    src,
    true);

  const unsigned int end_of_print_dst =
    dof_handler.n_dofs() > 9 ? 9 : dof_handler.n_dofs();

  deallog << "dst FEE: ";
  for (unsigned int i = 0; i < end_of_print_dst; ++i)
    deallog << dst[i] << ' ';
  deallog << std::endl;
}

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::hessians;
    deallog << "test_hessians_only" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<1>(i, dealii::FE_Q<1>(i), evaluation_flags, false);
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags, false);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags, false);
        test_hessians<1>(i, dealii::FE_DGQ<1>(i), evaluation_flags, false);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags, false);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags, false);
      }
    test_hessians<2>(3, dealii::FE_DGQ<2>(3), evaluation_flags, true);
  }
}
