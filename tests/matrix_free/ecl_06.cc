// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Compare the evaluation of a SIP Laplace operator with face- and element-
// centric loops on a hypercube, hypershell, and hyperball and p adaptivity.

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

template <int dim, typename Number, typename VectorizedArrayType>
void
test(const unsigned int geometry,
     const int          fe_degree,
     const bool         run_ecl,
     const unsigned int n_refinements = 1,
     const bool         print_vector  = true)
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

  tria.refine_global(n_refinements);

  FE_DGQ<dim>           fe1(fe_degree);
  FE_DGQ<dim>           fe2(fe_degree + 1);
  hp::FECollection<dim> fe(fe1, fe2);


  DoFHandler<dim> dof_handler(tria);
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        if (cell->center()[0] < 0.5)
          cell->set_active_fe_index(0);
        else
          cell->set_active_fe_index(1);
      }

  dof_handler.distribute_dofs(fe);

  MappingQ<dim>      mapping(1);
  hp::QCollection<1> quad;
  quad.push_back(QGauss<1>(fe_degree + 1));
  quad.push_back(QGauss<1>(fe_degree + 2));

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

  if (run_ecl)
    {
      AssertDimension(VectorizedArrayType::size(), 1);

      // note: catergorization and hp is not working together, i.e., we
      // cannot categorize the cells according to the boundary conditions
      // needed for ECL (TODO)

      additional_data.hold_all_faces_to_owned_cells        = true;
      additional_data.cell_vectorization_categories_strict = true;
      additional_data.mapping_update_flags_faces_by_cells =
        additional_data.mapping_update_flags_inner_faces |
        additional_data.mapping_update_flags_boundary_faces;

      // MatrixFreeTools::categorize_by_boundary_ids(tria, additional_data);
    }


  MF matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  VectorType src, dst;

  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  CosineFunction<dim, Number> function;
  VectorTools::interpolate(dof_handler, function, src);

  dst = 0.0;

  /**
   * Face-centric loop
   */
  matrix_free.template loop<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      const unsigned int active_fe_index =
        matrix_free.get_cell_iterator(range.first, 0)->active_fe_index();

      FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> phi(
        matrix_free, 0, 0, 0, active_fe_index, active_fe_index);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(EvaluationFlags::gradients);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_gradient(phi.get_gradient(q), q);
          phi.integrate(EvaluationFlags::gradients);
          phi.set_dof_values(dst);
        }
    },
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      const unsigned int active_fe_index_m =
        matrix_free.get_face_iterator(range.first, 0, true)
          .first->active_fe_index();
      const unsigned int active_fe_index_p =
        matrix_free.get_face_iterator(range.first, 0, false)
          .first->active_fe_index();
      const unsigned int active_q_index =
        std::max(active_fe_index_m, active_fe_index_p);

      FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> phi_m(
        matrix_free, true, 0, 0, 0, active_fe_index_m, active_q_index);
      FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> phi_p(
        matrix_free, false, 0, 0, 0, active_fe_index_p, active_q_index);

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
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      const unsigned int active_fe_index_m =
        matrix_free.get_face_iterator(range.first, 0, true)
          .first->active_fe_index();
      FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> phi_m(
        matrix_free, true, 0, 0, 0, active_fe_index_m, active_fe_index_m);
      for (unsigned int face = range.first; face < range.second; face++)
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
    dst.print(deallog.get_file_stream());

  if (run_ecl)
    {
      dst = 0.0;

      /**
       * Element-centric loop
       */
      matrix_free.template loop_cell_centric<VectorType, VectorType>(
        [&](const auto &, auto &dst, const auto &src, const auto range) {
          const unsigned int active_fe_index_m =
            matrix_free.get_cell_iterator(range.first, 0)->active_fe_index();

          FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> phi(
            matrix_free, 0, 0, 0, active_fe_index_m, active_fe_index_m);

          Table<2,
                std::shared_ptr<
                  FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>>>
            phi_ms(2, 2);

          for (unsigned int q = 0; q < quad.size(); ++q)
            for (unsigned int i = 0;
                 i < matrix_free.get_dof_handler().get_fe_collection().size();
                 ++i)
              phi_ms[i][q] = std::make_shared<
                FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>>(
                matrix_free, true, 0, 0, 0, i, q);

          Table<2,
                std::shared_ptr<
                  FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>>>
            phi_ps(2, 2);

          for (unsigned int q = 0; q < quad.size(); ++q)
            for (unsigned int i = 0;
                 i < matrix_free.get_dof_handler().get_fe_collection().size();
                 ++i)
              phi_ps[i][q] = std::make_shared<
                FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>>(
                matrix_free, false, 0, 0, 0, i, q);

          for (unsigned int cell = range.first; cell < range.second; ++cell)
            {
              phi.reinit(cell);
              phi.read_dof_values(src);
              phi.evaluate(EvaluationFlags::gradients);
              for (unsigned int q = 0; q < phi.n_q_points; ++q)
                phi.submit_gradient(phi.get_gradient(q), q);
              phi.integrate(EvaluationFlags::gradients);

              for (unsigned int face = 0;
                   face < GeometryInfo<dim>::faces_per_cell;
                   ++face)
                {
                  auto bids =
                    matrix_free.get_faces_by_cells_boundary_id(cell, face);

                  if (bids[0] != numbers::internal_face_boundary_id)
                    {
                      auto &phi_m =
                        *phi_ms[active_fe_index_m][active_fe_index_m];

                      phi_m.reinit(cell, face);
                      phi_m.read_dof_values(src);
                      phi_m.evaluate(EvaluationFlags::values |
                                     EvaluationFlags::gradients);
                      VectorizedArrayType sigmaF =
                        std::abs((phi_m.normal_vector(0) *
                                  phi_m.inverse_jacobian(0))[dim - 1]) *
                        Number(std::max(fe_degree, 1) * (fe_degree + 1.0)) *
                        2.0;

                      for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
                        {
                          VectorizedArrayType average_value =
                            phi_m.get_value(q);
                          VectorizedArrayType average_valgrad =
                            -phi_m.get_normal_derivative(q);
                          average_valgrad += average_value * sigmaF * 2.0;
                          phi_m.submit_normal_derivative(-average_value, q);
                          phi_m.submit_value(average_valgrad, q);
                        }

                      phi_m.integrate(EvaluationFlags::values |
                                      EvaluationFlags::gradients);

                      AssertDimension(phi.dofs_per_cell, phi_m.dofs_per_cell);

                      for (unsigned int q = 0; q < phi.dofs_per_cell; ++q)
                        phi.begin_dof_values()[q] +=
                          phi_m.begin_dof_values()[q];
                    }
                  else
                    {
                      const unsigned int active_fe_index_p =
                        matrix_free.get_cell_iterator(cell, 0)
                          ->neighbor(face)
                          ->active_fe_index();
                      const unsigned int active_q_index =
                        std::max(active_fe_index_m, active_fe_index_p);

                      auto &phi_m = *phi_ms[active_fe_index_m][active_q_index];
                      auto &phi_p = *phi_ps[active_fe_index_p][active_q_index];

                      phi_m.reinit(cell, face);
                      phi_p.reinit(cell, face);

                      phi_m.read_dof_values(src);
                      phi_m.evaluate(EvaluationFlags::values |
                                     EvaluationFlags::gradients);
                      phi_p.read_dof_values(src);
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

                      for (unsigned int q = 0; q < phi.dofs_per_cell; ++q)
                        phi.begin_dof_values()[q] +=
                          phi_m.begin_dof_values()[q];
                    }
                }

              phi.set_dof_values(dst);
            }
        },
        dst,
        src,
        false,
        MF::DataAccessOnFaces::gradients);

      deallog << dst.l2_norm() << std::endl;

      if (print_vector)
        dst.print(deallog.get_file_stream());

      std::cout << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  mpi_initlog();
  for (unsigned int i = 0; i < 3; ++i)
    {
      test<2, double, VectorizedArray<double>>(i, 2, false);
      test<2, double, VectorizedArray<double, 1>>(i, 2, true);
    }
}
