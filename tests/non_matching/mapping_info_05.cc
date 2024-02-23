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
do_flux_term_ecl(Integrator        &evaluator_m,
                 Integrator        &evaluator_p,
                 const Number2     &tau,
                 const unsigned int q)
{
  const auto gradient_m = evaluator_m.get_gradient(q);
  const auto gradient_p = evaluator_p.get_gradient(q);

  const auto value_m = evaluator_m.get_value(q);
  const auto value_p = evaluator_p.get_value(q);

  const auto normal = evaluator_m.normal_vector(q);

  const auto jump_value = (value_m - value_p) * normal;

  const auto central_flux_gradient = 0.5 * (gradient_m + gradient_p);

  const auto value_terms = normal * (central_flux_gradient - tau * jump_value);

  evaluator_m.submit_value(-value_terms, q);

  evaluator_m.submit_gradient(-0.5 * jump_value, q);
}

template <int dim>
void
test_dg_ecl(const unsigned int degree, const bool curved_mesh)
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
  additional_data.hold_all_faces_to_owned_cells = true;
  additional_data.mapping_update_flags_faces_by_cells =
    additional_data.mapping_update_flags_inner_faces |
    additional_data.mapping_update_flags_boundary_faces;

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

  matrix_free
    .template loop_cell_centric<LinearAlgebra::distributed::Vector<double>,
                                LinearAlgebra::distributed::Vector<double>>(
      [&](const auto &matrix_free,
          auto       &dst,
          const auto &src,
          const auto &range) {
        FEEvaluation<dim, -1>                  fe_eval(matrix_free);
        FEFaceEvaluation<dim, -1>              fe_eval_m(matrix_free, true);
        FEFaceEvaluation<dim, -1>              fe_eval_p(matrix_free, false);
        AlignedVector<VectorizedArray<double>> vec_solution_values_in_m(
          fe_eval.dofs_per_cell);
        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            fe_eval.reinit(cell);
            fe_eval.read_dof_values(src);
            for (unsigned int i = 0; i < fe_eval.dofs_per_cell; ++i)
              vec_solution_values_in_m[i] = fe_eval.begin_dof_values()[i];
            fe_eval.evaluate(EvaluationFlags::gradients);

            for (const unsigned int q : fe_eval.quadrature_point_indices())
              fe_eval.submit_gradient(fe_eval.get_gradient(q), q);

            fe_eval.integrate(EvaluationFlags::gradients);

            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
              {
                // ask for boundary ids of face
                const auto boundary_ids =
                  matrix_free.get_faces_by_cells_boundary_id(cell, f);

                // only internal faces have a neighbor, setup a mask
                std::bitset<n_lanes>    mask;
                VectorizedArray<double> fluxes = 0.;
                for (unsigned int v = 0; v < n_lanes; ++v)
                  {
                    mask[v] =
                      boundary_ids[v] == numbers::internal_face_boundary_id;
                    fluxes[v] = mask[v] == true ? 1. : 0.;
                  }

                fe_eval_m.reinit(cell, f);
                fe_eval_p.reinit(cell, f);

                fe_eval_p.read_dof_values(src, 0, mask);

                fe_eval_m.evaluate(vec_solution_values_in_m.data(),
                                   EvaluationFlags::values |
                                     EvaluationFlags::gradients);
                fe_eval_p.evaluate(EvaluationFlags::values |
                                   EvaluationFlags::gradients);

                for (const auto q : fe_eval_m.quadrature_point_indices())
                  {
                    do_flux_term_ecl(fe_eval_m, fe_eval_p, 1.0, q);

                    // clear lanes where face at boundary
                    fe_eval_m.begin_values()[q] *= fluxes;
                    for (unsigned int d = 0; d < dim; ++d)
                      fe_eval_m.begin_gradients()[q * dim + d] *= fluxes;
                  }

                fe_eval_m.integrate(EvaluationFlags::values |
                                      EvaluationFlags::gradients,
                                    fe_eval.begin_dof_values(),
                                    true);
              }

            fe_eval.distribute_local_to_global(dst);
          }
      },
      dst,
      src,
      true);

  QGauss<dim>                  quad_cell(n_q_points);
  QGauss<dim - 1>              quad_face(n_q_points);
  std::vector<Quadrature<dim>> quad_vec_cells;
  quad_vec_cells.reserve(
    (matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches()) *
    n_lanes);
  std::vector<std::vector<Quadrature<dim - 1>>> quad_vec_faces(
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

        for (const auto f : GeometryInfo<dim>::face_indices())
          {
            (void)f;
            quad_vec_faces[cell_batch * n_lanes + v].push_back(quad_face);
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
  mapping_info_faces.reinit_faces(vector_accessors, quad_vec_faces);

  FEPointEvaluation<1, dim, dim, double>     fe_peval(mapping_info_cells, fe);
  FEFacePointEvaluation<1, dim, dim, double> fe_peval_m(mapping_info_faces, fe);
  FEFacePointEvaluation<1, dim, dim, double> fe_peval_p(mapping_info_faces, fe);

  matrix_free
    .template loop_cell_centric<LinearAlgebra::distributed::Vector<double>,
                                LinearAlgebra::distributed::Vector<double>>(
      [&](const auto &matrix_free,
          auto       &dst,
          const auto &src,
          const auto &range) {
        FEEvaluation<dim, -1>                  fe_eval(matrix_free);
        FEFaceEvaluation<dim, -1>              fe_eval_m(matrix_free, true);
        FEFaceEvaluation<dim, -1>              fe_eval_p(matrix_free, false);
        AlignedVector<VectorizedArray<double>> vec_solution_values_in_m(
          fe_eval.dofs_per_cell);
        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            fe_eval.reinit(cell);
            fe_eval.read_dof_values(src);

            for (unsigned int i = 0; i < fe_eval.dofs_per_cell; ++i)
              vec_solution_values_in_m[i] = fe_eval.begin_dof_values()[i];

            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                fe_peval.reinit(cell * n_lanes + v);
                fe_peval.evaluate(StridedArrayView<const double, n_lanes>(
                                    &vec_solution_values_in_m[0][v],
                                    fe.dofs_per_cell),
                                  EvaluationFlags::gradients);
                for (const unsigned int q : fe_peval.quadrature_point_indices())
                  fe_peval.submit_gradient(fe_peval.get_gradient(q), q);
                fe_peval.integrate(StridedArrayView<double, n_lanes>(
                                     &fe_eval.begin_dof_values()[0][v],
                                     fe.dofs_per_cell),
                                   EvaluationFlags::gradients);
              }

            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
              {
                fe_eval_m.reinit(cell, f);

                fe_eval_m.project_to_face(&vec_solution_values_in_m.begin()[0],
                                          EvaluationFlags::values |
                                            EvaluationFlags::gradients);

                // ask for boundary ids of face
                const auto boundary_ids =
                  matrix_free.get_faces_by_cells_boundary_id(cell, f);

                // only internal faces have a neighbor, setup a mask
                std::bitset<n_lanes> mask;
                for (unsigned int v = 0; v < n_lanes; ++v)
                  mask[v] =
                    boundary_ids[v] == numbers::internal_face_boundary_id;

                fe_eval_p.reinit(cell, f);
                fe_eval_p.read_dof_values(src, 0, mask);

                fe_eval_p.project_to_face(EvaluationFlags::values |
                                          EvaluationFlags::gradients);

                const auto &cell_indices_p = fe_eval_p.get_cell_ids();

                for (unsigned int v = 0; v < n_lanes; ++v)
                  {
                    if (mask[v] == false)
                      {
                        for (unsigned int i = 0;
                             i < 2 * fe_eval_m.get_dofs_projected_to_face();
                             ++i)
                          fe_eval_m.get_scratch_data().begin()[i][v] = 0.;
                        continue;
                      }

                    fe_peval_m.reinit(cell * n_lanes + v, f);

                    fe_peval_p.reinit(cell_indices_p[v],
                                      fe_eval_p.get_face_no(v));

                    fe_peval_m.evaluate_in_face(
                      &fe_eval_m.get_scratch_data().begin()[0][v],
                      EvaluationFlags::values | EvaluationFlags::gradients);
                    fe_peval_p.evaluate_in_face(
                      &fe_eval_p.get_scratch_data().begin()[0][v],
                      EvaluationFlags::values | EvaluationFlags::gradients);

                    for (const unsigned int q :
                         fe_peval_m.quadrature_point_indices())
                      do_flux_term_ecl(fe_peval_m, fe_peval_p, 1.0, q);

                    fe_peval_m.integrate_in_face(
                      &fe_eval_m.get_scratch_data().begin()[0][v],
                      EvaluationFlags::values | EvaluationFlags::gradients);
                  }

                fe_eval_m.collect_from_face(EvaluationFlags::values |
                                              EvaluationFlags::gradients,
                                            fe_eval.begin_dof_values(),
                                            true);
              }

            fe_eval.distribute_local_to_global(dst);
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

  test_dg_ecl<2>(1, false);
  deallog << std::endl;
  test_dg_ecl<2>(2, false);
  deallog << std::endl;
  test_dg_ecl<3>(1, false);
  deallog << std::endl;
  test_dg_ecl<3>(2, false);
  deallog << std::endl;
  test_dg_ecl<2>(1, true);
  deallog << std::endl;
  test_dg_ecl<2>(2, true);
  // TODO: fix face orientation
  // deallog << std::endl;
  // test_dg_ecl<3>(1, true);
  // deallog << std::endl;
  // test_dg_ecl<3>(2, true);
}
