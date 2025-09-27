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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// Adapted from ecl.h

template <int dim, typename Number>
class CosineFunction : public Function<dim, Number>
{
public:
  Number
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;
    return std::cos(2 * numbers::PI * p[0]);
  }
};

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const unsigned int geometry      = 0,
     const unsigned int n_refinements = 1,
     const bool         print_vector  = true,
     const MPI_Comm     comm          = MPI_COMM_SELF)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  if (geometry == 0)
    GridGenerator::hyper_cube(tria, 0.0, 1.0, true);
  else if (geometry == 1)
    GridGenerator::hyper_shell(tria, Point<dim>(), 0.5, 1.0);
  else if (geometry == 2)
    GridGenerator::hyper_ball(tria);
  else
    DEAL_II_NOT_IMPLEMENTED();

  tria.reset_all_manifolds();

  if (false)
    {
      std::vector<dealii::GridTools::PeriodicFacePair<
        typename dealii::Triangulation<dim>::cell_iterator>>
        periodic_faces;

      if (dim >= 1)
        dealii::GridTools::collect_periodic_faces(
          tria, 0, 1, 0, periodic_faces);

      if (dim >= 2)
        dealii::GridTools::collect_periodic_faces(
          tria, 2, 3, 1, periodic_faces);

      if (dim >= 3)
        dealii::GridTools::collect_periodic_faces(
          tria, 4, 5, 2, periodic_faces);

      tria.add_periodicity(periodic_faces);
    }

  tria.refine_global(n_refinements);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(n_points);

  AffineConstraints<Number> constraint;

  using MF = MatrixFree<dim, Number, VectorizedArrayType>;

  typename MF::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients;
  additional_data.mapping_update_flags_boundary_faces =
    update_values | update_gradients;
  additional_data.tasks_parallel_scheme =
    MF::AdditionalData::TasksParallelScheme::none;
  additional_data.mapping_update_flags_faces_by_cells =
    update_values | update_gradients;
  additional_data.hold_all_faces_to_owned_cells = true;
  additional_data.communicator_sm               = comm;

  MF matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  VectorType src, dst;

  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  FEEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType> phi(
    matrix_free);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_m(matrix_free, true);
  FEFaceEvaluation<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_p(matrix_free, false);


  CosineFunction<dim, Number> function;
  VectorTools::interpolate(dof_handler, function, src);

  dst = 0.0;

  /**
   * Face-centric loop
   */
  matrix_free.template loop<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(EvaluationFlags::gradients);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(EvaluationFlags::gradients);
          phi.distribute_local_to_global(dst);
        }
    },
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          phi_m.reinit(face);
          phi_p.reinit(face);

          phi_m.read_dof_values(src);
          phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
          phi_p.read_dof_values(src);
          phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
          VectorizedArrayType sigmaF =
            (std::abs(
               (phi_m.normal_vector(0) * phi_m.inverse_jacobian(0))[dim - 1]) +
             std::abs(
               (phi_m.normal_vector(0) * phi_p.inverse_jacobian(0))[dim - 1])) *
            (Number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

          for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
            {
              VectorizedArrayType average_value =
                (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
              VectorizedArrayType average_valgrad =
                phi_m.get_normal_derivative(q) + phi_p.get_normal_derivative(q);
              average_valgrad =
                average_value * 2. * sigmaF - average_valgrad * 0.5;
              phi_m.submit_normal_derivative(-average_value, q);
              phi_p.submit_normal_derivative(-average_value, q);
              phi_m.submit_value(average_valgrad, q);
              phi_p.submit_value(-average_valgrad, q);
            }
          phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
          phi_m.distribute_local_to_global(dst);
          phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
          phi_p.distribute_local_to_global(dst);
        }
    },
    [&](const auto &, auto &dst, const auto &src, const auto face_range) {
      for (unsigned int face = face_range.first; face < face_range.second;
           face++)
        {
          phi_m.reinit(face);
          phi_m.read_dof_values(src);
          phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
          VectorizedArrayType sigmaF =
            std::abs(
              (phi_m.normal_vector(0) * phi_m.inverse_jacobian(0))[dim - 1]) *
            Number(std::max(fe_degree, 1) * (fe_degree + 1.0)) * 2.0;

          for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
            {
              VectorizedArrayType average_value = phi_m.get_value(q);
              VectorizedArrayType average_valgrad =
                -phi_m.get_normal_derivative(q);
              average_valgrad += average_value * sigmaF * 2.0;
              phi_m.submit_normal_derivative(-average_value, q);
              phi_m.submit_value(average_valgrad, q);
            }

          phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
          phi_m.distribute_local_to_global(dst);
        }
    },
    dst,
    src,
    false,
    MF::DataAccessOnFaces::gradients,
    MF::DataAccessOnFaces::gradients);

  deallog << dst.l2_norm() << std::endl;

  if (print_vector)
    dst.print(deallog.get_file_stream(), 10);

  dst = 0.0;

  /**
   * Element-centric loop
   */
  matrix_free.template loop_cell_centric<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(EvaluationFlags::gradients);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(EvaluationFlags::gradients);

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            {
              const auto bids =
                matrix_free.get_faces_by_cells_boundary_id(cell, face);

              std::bitset<VectorizedArrayType::size()> mask;
              for (unsigned int v = 0;
                   v < matrix_free.n_active_entries_per_cell_batch(cell);
                   ++v)
                if (bids[v] == numbers::internal_face_boundary_id)
                  mask[v] = true;

              // boundary
              {
                phi_m.reinit(cell, face);
                phi_m.read_dof_values(src);
                phi_m.evaluate(EvaluationFlags::values |
                               EvaluationFlags::gradients);
                VectorizedArrayType sigmaF =
                  std::abs((phi_m.normal_vector(0) *
                            phi_m.inverse_jacobian(0))[dim - 1]) *
                  Number(std::max(fe_degree, 1) * (fe_degree + 1.0)) * 2.0;

                for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                  {
                    VectorizedArrayType average_value = phi_m.get_value(q);
                    VectorizedArrayType average_valgrad =
                      -phi_m.get_normal_derivative(q);
                    average_valgrad += average_value * sigmaF * 2.0;
                    phi_m.submit_normal_derivative(-average_value, q);
                    phi_m.submit_value(average_valgrad, q);
                  }

                phi_m.integrate(EvaluationFlags::values |
                                EvaluationFlags::gradients);

                for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                  for (unsigned int q = 0; q < phi.dofs_per_cell; ++q)
                    phi.begin_dof_values()[q][v] +=
                      mask[v] == false ? phi_m.begin_dof_values()[q][v] : 0.;
              }

              // interior face
              {
                phi_m.reinit(cell, face);
                phi_p.reinit(cell, face);

                phi_m.read_dof_values(src);
                phi_m.evaluate(EvaluationFlags::values |
                               EvaluationFlags::gradients);

                phi_p.read_dof_values(src, 0, mask);
                phi_p.evaluate(EvaluationFlags::values |
                               EvaluationFlags::gradients);

                VectorizedArrayType sigmaF =
                  (std::abs((phi_m.normal_vector(0) *
                             phi_m.inverse_jacobian(0))[dim - 1]) +
                   std::abs((phi_m.normal_vector(0) *
                             phi_p.inverse_jacobian(0))[dim - 1])) *
                  (Number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

                for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                  {
                    VectorizedArrayType average_value =
                      (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
                    VectorizedArrayType average_valgrad =
                      phi_m.get_normal_derivative(q) +
                      phi_p.get_normal_derivative(q);
                    average_valgrad =
                      average_value * 2. * sigmaF - average_valgrad * 0.5;
                    phi_m.submit_normal_derivative(-average_value, q);
                    phi_m.submit_value(average_valgrad, q);
                  }
                phi_m.integrate(EvaluationFlags::values |
                                EvaluationFlags::gradients);

                for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                  for (unsigned int q = 0; q < phi.dofs_per_cell; ++q)
                    phi.begin_dof_values()[q][v] +=
                      mask[v] == true ? phi_m.begin_dof_values()[q][v] : 0.;
              }
            }

          phi.distribute_local_to_global(dst);
        }
    },
    dst,
    src,
    false,
    MF::DataAccessOnFaces::gradients);

  deallog << dst.l2_norm() << std::endl;

  if (print_vector)
    dst.print(deallog.get_file_stream(), 10);
}
