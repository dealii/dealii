// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Test the NonMatching::MappingInfo class together with FEPointEvaluation and
 * compare to FEEvaluation
 */

#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <typename Integrator, typename Number2>
void
do_flux_term(Integrator        &evaluator_m,
             Integrator        &evaluator_p,
             const Number2     &tau,
             const unsigned int q)
{
  const auto normal_gradient_m = evaluator_m.get_normal_derivative(q);
  const auto normal_gradient_p = evaluator_p.get_normal_derivative(q);

  const auto value_m = evaluator_m.get_value(q);
  const auto value_p = evaluator_p.get_value(q);

  const auto jump_value = value_m - value_p;

  const auto central_flux_gradient =
    0.5 * (normal_gradient_m + normal_gradient_p);

  const auto value_terms = central_flux_gradient - tau * jump_value;

  evaluator_m.submit_value(-value_terms, q);
  evaluator_p.submit_value(value_terms, q);

  const auto gradient_terms = -0.5 * jump_value;

  evaluator_m.submit_normal_derivative(gradient_terms, q);
  evaluator_p.submit_normal_derivative(gradient_terms, q);
}

template <int dim>
void
test_dg_fcl(const unsigned int degree, const bool curved_mesh)
{
  constexpr unsigned int n_lanes = VectorizedArray<double>::size();

  const unsigned int n_q_points = degree + 1;

  Triangulation<dim> tria;

  if (curved_mesh && dim > 1)
    GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1, 6);
  else
    GridGenerator::subdivided_hyper_cube(tria, 2, 0, 1);

  tria.refine_global(1);

  FE_DGQ<dim>     fe(degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  MappingQGeneric<dim> mapping(degree);

  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  typename MatrixFree<dim>::AdditionalData additional_data;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients;
  additional_data.mapping_update_flags_boundary_faces =
    update_values | update_gradients;

  MatrixFree<dim> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, QGauss<1>(n_q_points), additional_data);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Working with " << fe.get_name() << " and "
            << dof_handler.n_dofs() << " dofs" << std::endl;

  LinearAlgebra::distributed::Vector<double> src, dst, dst2;
  matrix_free.initialize_dof_vector(src);
  for (auto &v : src)
    v = static_cast<double>(rand()) / RAND_MAX;

  matrix_free.initialize_dof_vector(dst);
  matrix_free.initialize_dof_vector(dst2);

  matrix_free.template loop<LinearAlgebra::distributed::Vector<double>,
                            LinearAlgebra::distributed::Vector<double>>(
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEEvaluation<dim, -1> fe_eval(matrix_free);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          fe_eval.reinit(cell);
          fe_eval.gather_evaluate(src, EvaluationFlags::gradients);
          for (const unsigned int q : fe_eval.quadrature_point_indices())
            fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
          fe_eval.integrate_scatter(EvaluationFlags::gradients, dst);
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1> fe_eval_m(matrix_free, true);
      FEFaceEvaluation<dim, -1> fe_eval_p(matrix_free, false);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          fe_eval_m.reinit(face);
          fe_eval_p.reinit(face);
          fe_eval_m.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
          fe_eval_p.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
          for (unsigned int q = 0; q < fe_eval_m.n_q_points; ++q)
            do_flux_term(fe_eval_m, fe_eval_p, 1.0, q);
          fe_eval_m.integrate_scatter(EvaluationFlags::values |
                                        EvaluationFlags::gradients,
                                      dst);
          fe_eval_p.integrate_scatter(EvaluationFlags::values |
                                        EvaluationFlags::gradients,
                                      dst);
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1> fe_eval_m(matrix_free, true);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          fe_eval_m.reinit(face);
          fe_eval_m.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
          for (unsigned int q = 0; q < fe_eval_m.n_q_points; ++q)
            {
              const auto value    = fe_eval_m.get_value(q);
              const auto gradient = fe_eval_m.get_gradient(q);

              fe_eval_m.submit_value(gradient * fe_eval_m.normal_vector(q) +
                                       value,
                                     q);
              fe_eval_m.submit_gradient(value * fe_eval_m.normal_vector(q), q);
            }
          fe_eval_m.integrate_scatter(EvaluationFlags::values |
                                        EvaluationFlags::gradients,
                                      dst);
        }
    },
    dst,
    src,
    true);

  QGauss<dim>                  quad_cell(n_q_points);
  std::vector<Quadrature<dim>> quad_vec_cells;
  quad_vec_cells.reserve(
    (matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches()) *
    n_lanes);


  std::vector<typename DoFHandler<dim>::cell_iterator> vector_accessors;
  vector_accessors.reserve(
    (matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches()) *
    n_lanes);
  for (unsigned int cell_batch = 0;
       cell_batch <
       matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
       ++cell_batch)
    for (unsigned int v = 0; v < n_lanes; ++v)
      {
        if (v < matrix_free.n_active_entries_per_cell_batch(cell_batch))
          vector_accessors.push_back(
            matrix_free.get_cell_iterator(cell_batch, v));
        else
          vector_accessors.push_back(
            matrix_free.get_cell_iterator(cell_batch, 0));

        quad_vec_cells.push_back(quad_cell);
      }

  QGauss<dim - 1>                  quad_face(n_q_points);
  std::vector<Quadrature<dim - 1>> quad_vec_faces;
  quad_vec_faces.reserve((matrix_free.n_inner_face_batches() +
                          matrix_free.n_boundary_face_batches()) *
                         n_lanes);
  std::vector<std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>>
    vector_face_accessors_m;
  vector_face_accessors_m.reserve((matrix_free.n_inner_face_batches() +
                                   matrix_free.n_boundary_face_batches()) *
                                  n_lanes);
  // fill container for inner face batches
  unsigned int face_batch = 0;
  for (; face_batch < matrix_free.n_inner_face_batches(); ++face_batch)
    {
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
            vector_face_accessors_m.push_back(
              matrix_free.get_face_iterator(face_batch, v));
          else
            vector_face_accessors_m.push_back(
              matrix_free.get_face_iterator(face_batch, 0));

          quad_vec_faces.push_back(quad_face);
        }
    }
  // and boundary face batches
  for (; face_batch < (matrix_free.n_inner_face_batches() +
                       matrix_free.n_boundary_face_batches());
       ++face_batch)
    {
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
            vector_face_accessors_m.push_back(
              matrix_free.get_face_iterator(face_batch, v));
          else
            vector_face_accessors_m.push_back(
              matrix_free.get_face_iterator(face_batch, 0));

          quad_vec_faces.push_back(quad_face);
        }
    }

  NonMatching::MappingInfo<dim> mapping_info_cells(mapping,
                                                   update_gradients |
                                                     update_JxW_values);
  NonMatching::MappingInfo<dim> mapping_info_faces(mapping,
                                                   update_values |
                                                     update_gradients |
                                                     update_JxW_values |
                                                     update_normal_vectors);

  mapping_info_cells.reinit_cells(vector_accessors, quad_vec_cells);
  mapping_info_faces.reinit_faces(vector_face_accessors_m, quad_vec_faces);

  FEPointEvaluation<1, dim, dim, double>     fe_peval(mapping_info_cells, fe);
  FEFacePointEvaluation<1, dim, dim, double> fe_peval_m(mapping_info_faces,
                                                        fe,
                                                        true);
  FEFacePointEvaluation<1, dim, dim, double> fe_peval_p(mapping_info_faces,
                                                        fe,
                                                        false);

  matrix_free.template loop<LinearAlgebra::distributed::Vector<double>,
                            LinearAlgebra::distributed::Vector<double>>(
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEEvaluation<dim, -1> fe_eval(matrix_free);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          fe_eval.reinit(cell);
          fe_eval.read_dof_values(src);
          for (unsigned int v = 0; v < n_lanes; ++v)
            {
              fe_peval.reinit(cell * n_lanes + v);
              fe_peval.evaluate(StridedArrayView<const double, n_lanes>(
                                  &fe_eval.begin_dof_values()[0][v],
                                  fe.dofs_per_cell),
                                EvaluationFlags::gradients);
              for (const unsigned int q : fe_peval.quadrature_point_indices())
                fe_peval.submit_gradient(fe_peval.get_gradient(q), q);
              fe_peval.integrate(StridedArrayView<double, n_lanes>(
                                   &fe_eval.begin_dof_values()[0][v],
                                   fe.dofs_per_cell),
                                 EvaluationFlags::gradients);
            }
          fe_eval.distribute_local_to_global(dst);
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1> fe_eval_m(matrix_free, true);
      FEFaceEvaluation<dim, -1> fe_eval_p(matrix_free, false);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          fe_eval_m.reinit(face);
          fe_eval_p.reinit(face);

          fe_eval_m.read_dof_values(src);
          fe_eval_p.read_dof_values(src);

          fe_eval_m.project_to_face(EvaluationFlags::values |
                                    EvaluationFlags::gradients);
          fe_eval_p.project_to_face(EvaluationFlags::values |
                                    EvaluationFlags::gradients);

          for (unsigned int v = 0; v < n_lanes; ++v)
            {
              fe_peval_m.reinit(face * n_lanes + v);
              fe_peval_p.reinit(face * n_lanes + v);
              fe_peval_m.evaluate_in_face(
                &fe_eval_m.get_scratch_data().begin()[0][v],
                EvaluationFlags::values | EvaluationFlags::gradients);
              fe_peval_p.evaluate_in_face(
                &fe_eval_p.get_scratch_data().begin()[0][v],
                EvaluationFlags::values | EvaluationFlags::gradients);
              for (const unsigned int q : fe_peval_m.quadrature_point_indices())
                do_flux_term(fe_peval_m, fe_peval_p, 1.0, q);
              fe_peval_m.integrate_in_face(
                &fe_eval_m.get_scratch_data().begin()[0][v],
                EvaluationFlags::values | EvaluationFlags::gradients);
              fe_peval_p.integrate_in_face(
                &fe_eval_p.get_scratch_data().begin()[0][v],
                EvaluationFlags::values | EvaluationFlags::gradients);
            }

          fe_eval_m.collect_from_face(EvaluationFlags::values |
                                      EvaluationFlags::gradients);
          fe_eval_p.collect_from_face(EvaluationFlags::values |
                                      EvaluationFlags::gradients);

          fe_eval_m.distribute_local_to_global(dst);
          fe_eval_p.distribute_local_to_global(dst);
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1> fe_eval_m(matrix_free, true);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          fe_eval_m.reinit(face);

          fe_eval_m.read_dof_values(src);

          fe_eval_m.project_to_face(EvaluationFlags::values |
                                    EvaluationFlags::gradients);

          for (unsigned int v = 0; v < n_lanes; ++v)
            {
              fe_peval_m.reinit(face * n_lanes + v);
              fe_peval_m.evaluate_in_face(
                &fe_eval_m.get_scratch_data().begin()[0][v],
                EvaluationFlags::values | EvaluationFlags::gradients);
              for (const unsigned int q : fe_peval_m.quadrature_point_indices())
                {
                  const auto value    = fe_peval_m.get_value(q);
                  const auto gradient = fe_peval_m.get_gradient(q);

                  fe_peval_m.submit_value(
                    gradient * fe_peval_m.normal_vector(q) + value, q);
                  fe_peval_m.submit_gradient(value *
                                               fe_peval_m.normal_vector(q),
                                             q);
                }
              fe_peval_m.integrate_in_face(
                &fe_eval_m.get_scratch_data().begin()[0][v],
                EvaluationFlags::values | EvaluationFlags::gradients);
            }

          fe_eval_m.collect_from_face(EvaluationFlags::values |
                                      EvaluationFlags::gradients);

          fe_eval_m.distribute_local_to_global(dst);
        }
    },
    dst2,
    src,
    true);


  dst2 -= dst;
  const double error = dst2.l2_norm() / dst.l2_norm();
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "FEPointEvaluation verification: " << error << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();

  test_dg_fcl<2>(1, false);
  deallog << std::endl;
  test_dg_fcl<2>(2, false);
  deallog << std::endl;
  test_dg_fcl<3>(1, false);
  deallog << std::endl;
  test_dg_fcl<3>(2, false);
  deallog << std::endl;
  test_dg_fcl<2>(1, true);
  deallog << std::endl;
  test_dg_fcl<2>(2, true);
  // TODO: fix face orientation
  // deallog << std::endl;
  // test_dg_fcl<3>(1, true);
  // deallog << std::endl;
  // test_dg_fcl<3>(2, true);
}
