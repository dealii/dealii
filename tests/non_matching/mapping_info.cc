// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
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
 * compare to NonMatching::FEValues
 */

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mapping_info.h>
#include <deal.II/non_matching/mesh_classifier.h>
#include <deal.II/non_matching/quadrature_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test(const bool filtered_compression)
{
  constexpr unsigned int degree = 1;

  FE_Q<dim> fe_q(degree);

  Triangulation<dim> tria;

  MappingQ<dim> mapping(degree);

  GridGenerator::subdivided_hyper_cube(tria, 4);

  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe_q);

  Functions::SignedDistance::Sphere<dim> level_set;

  Vector<double> level_set_vec(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, level_set, level_set_vec);

  NonMatching::MeshClassifier<dim> mesh_classifier(dof_handler, level_set_vec);
  mesh_classifier.reclassify();

  hp::QCollection<1> q_collection((QGauss<1>(degree)));

  NonMatching::DiscreteQuadratureGenerator<dim> quadrature_generator(
    q_collection, dof_handler, level_set_vec);

  NonMatching::DiscreteFaceQuadratureGenerator<dim> face_quadrature_generator(
    q_collection, dof_handler, level_set_vec);

  // FEPointEvaluation
  NonMatching::MappingInfo<dim> mapping_info_cell(
    mapping, update_values | update_gradients | update_JxW_values);

  NonMatching::MappingInfo<dim> mapping_info_surface(mapping,
                                                     update_values |
                                                       update_gradients |
                                                       update_JxW_values |
                                                       update_normal_vectors);

  NonMatching::MappingInfo<dim> mapping_info_faces(mapping,
                                                   update_values |
                                                     update_gradients |
                                                     update_JxW_values |
                                                     update_normal_vectors);

  auto physical =
    [&](const typename Triangulation<dim>::active_cell_iterator &i) {
      return (mesh_classifier.location_to_level_set(i) ==
              NonMatching::LocationToLevelSet::intersected) ||
             (mesh_classifier.location_to_level_set(i) ==
              NonMatching::LocationToLevelSet::inside);
    };

  if (filtered_compression)
    {
      const auto physical_active_cell_iterators =
        tria.active_cell_iterators() | physical;

      unsigned int n_physical_cells = 0;
      for (const auto &cell : physical_active_cell_iterators)
        {
          ++n_physical_cells;
        }

      std::vector<Quadrature<dim>> quad_vec_cell;
      quad_vec_cell.reserve(n_physical_cells);
      std::vector<NonMatching::ImmersedSurfaceQuadrature<dim>> quad_vec_surface;
      quad_vec_surface.reserve(n_physical_cells);
      std::vector<std::vector<Quadrature<dim - 1>>> quad_vec_faces(
        n_physical_cells);

      unsigned int counter = 0;
      for (const auto &cell : physical_active_cell_iterators)
        {
          quadrature_generator.generate(cell);
          quad_vec_cell.push_back(quadrature_generator.get_inside_quadrature());
          quad_vec_surface.push_back(
            quadrature_generator.get_surface_quadrature());

          for (auto f : cell->face_indices())
            {
              face_quadrature_generator.generate(cell, f);
              quad_vec_faces[counter].push_back(
                face_quadrature_generator.get_inside_quadrature());
            }

          ++counter;
        }

      mapping_info_cell.reinit_cells(physical_active_cell_iterators,
                                     quad_vec_cell,
                                     tria.n_active_cells());
      mapping_info_surface.reinit_surface(physical_active_cell_iterators,
                                          quad_vec_surface,
                                          tria.n_active_cells());
      mapping_info_faces.reinit_faces(physical_active_cell_iterators,
                                      quad_vec_faces,
                                      tria.n_active_cells());
    }
  else
    {
      std::vector<Quadrature<dim>>                             quad_vec_cell;
      std::vector<NonMatching::ImmersedSurfaceQuadrature<dim>> quad_vec_surface;
      std::vector<std::vector<Quadrature<dim - 1>>>            quad_vec_faces(
        tria.n_active_cells());
      for (const auto &cell : tria.active_cell_iterators())
        {
          quadrature_generator.generate(cell);
          quad_vec_cell.push_back(quadrature_generator.get_inside_quadrature());
          quad_vec_surface.push_back(
            quadrature_generator.get_surface_quadrature());

          for (auto f : cell->face_indices())
            {
              face_quadrature_generator.generate(cell, f);
              quad_vec_faces[cell->active_cell_index()].push_back(
                face_quadrature_generator.get_inside_quadrature());
            }
        }

      mapping_info_cell.reinit_cells(tria.active_cell_iterators(),
                                     quad_vec_cell);
      mapping_info_surface.reinit_surface(tria.active_cell_iterators(),
                                          quad_vec_surface);
      mapping_info_faces.reinit_faces(tria.active_cell_iterators(),
                                      quad_vec_faces);
    }

  Vector<double> src(dof_handler.n_dofs()), dst_cell(dof_handler.n_dofs()),
    dst_surface(dof_handler.n_dofs()), dst_faces(dof_handler.n_dofs());

  for (auto &v : src)
    v = random_value<double>();

  FEPointEvaluation<1, dim, dim, double> fe_point_cell(mapping_info_cell, fe_q);
  FEPointEvaluation<1, dim, dim, double> fe_point_surface(mapping_info_surface,
                                                          fe_q);
  FEFacePointEvaluation<1, dim, dim, double> fe_point_faces_m(
    mapping_info_faces, fe_q);
  FEFacePointEvaluation<1, dim, dim, double> fe_point_faces_p(
    mapping_info_faces, fe_q);

  std::vector<double> solution_values_in(fe_q.dofs_per_cell);
  std::vector<double> solution_values_neighbor_in(fe_q.dofs_per_cell);
  std::vector<double> solution_values_cell_out(fe_q.dofs_per_cell);
  std::vector<double> solution_values_surface_out(fe_q.dofs_per_cell);
  std::vector<double> solution_values_faces_out(fe_q.dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!(physical(cell)))
        continue;

      cell->get_dof_values(src,
                           solution_values_in.begin(),
                           solution_values_in.end());

      auto test_fe_point = [&](FEPointEvaluation<1, dim, dim, double> &fe_point,
                               std::vector<double> &solution_values_out) {
        fe_point.reinit(cell->active_cell_index());

        fe_point.evaluate(solution_values_in,
                          EvaluationFlags::values | EvaluationFlags::gradients);

        for (const auto q : fe_point.quadrature_point_indices())
          {
            fe_point.submit_value(fe_point.get_value(q), q);
            fe_point.submit_gradient(fe_point.get_gradient(q), q);
          }

        fe_point.integrate(solution_values_out,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);
      };

      test_fe_point(fe_point_cell, solution_values_cell_out);
      test_fe_point(fe_point_surface, solution_values_surface_out);

      for (const auto f : cell->face_indices())
        {
          if (cell->at_boundary(f) || !physical(cell->neighbor(f)))
            continue;

          fe_point_faces_m.reinit(cell->active_cell_index(), f);
          fe_point_faces_p.reinit(cell->neighbor(f)->active_cell_index(),
                                  cell->neighbor_of_neighbor(f));

          cell->neighbor(f)->get_dof_values(src,
                                            solution_values_neighbor_in.begin(),
                                            solution_values_neighbor_in.end());

          fe_point_faces_m.evaluate(solution_values_in,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);

          fe_point_faces_p.evaluate(solution_values_neighbor_in,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);

          for (const auto q : fe_point_faces_m.quadrature_point_indices())
            {
              fe_point_faces_m.submit_value(fe_point_faces_m.get_value(q) -
                                              fe_point_faces_p.get_value(q),
                                            q);
              fe_point_faces_m.submit_gradient(
                fe_point_faces_m.get_gradient(q) -
                  fe_point_faces_p.get_gradient(q),
                q);
            }

          fe_point_faces_m.integrate(solution_values_faces_out,
                                     EvaluationFlags::values |
                                       EvaluationFlags::gradients);

          cell->distribute_local_to_global(
            Vector<double>(solution_values_faces_out.begin(),
                           solution_values_faces_out.end()),
            dst_faces);
        }

      cell->distribute_local_to_global(
        Vector<double>(solution_values_cell_out.begin(),
                       solution_values_cell_out.end()),
        dst_cell);

      cell->distribute_local_to_global(
        Vector<double>(solution_values_surface_out.begin(),
                       solution_values_surface_out.end()),
        dst_surface);
    }


  // FEValues
  const QGauss<1> quadrature_1D(degree);

  NonMatching::RegionUpdateFlags region_update_flags;
  region_update_flags.inside = update_values | update_gradients |
                               update_JxW_values | update_quadrature_points;
  region_update_flags.surface = update_values | update_gradients |
                                update_JxW_values | update_quadrature_points |
                                update_normal_vectors;

  hp::FECollection<dim> fe_collection(fe_q);

  NonMatching::FEValues<dim> non_matching_fe_values(fe_collection,
                                                    quadrature_1D,
                                                    region_update_flags,
                                                    mesh_classifier,
                                                    dof_handler,
                                                    level_set_vec);

  NonMatching::FEInterfaceValues<dim> non_matching_fe_interface_values(
    fe_collection,
    quadrature_1D,
    region_update_flags,
    mesh_classifier,
    dof_handler,
    level_set_vec);

  Vector<double> dst_cell_2(dof_handler.n_dofs()),
    dst_surface_2(dof_handler.n_dofs()), dst_faces_2(dof_handler.n_dofs());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      non_matching_fe_values.reinit(cell);

      cell->get_dof_values(src,
                           solution_values_in.begin(),
                           solution_values_in.end());

      const auto &inside_fe_values =
        non_matching_fe_values.get_inside_fe_values();

      const auto &surface_fe_values =
        non_matching_fe_values.get_surface_fe_values();


      auto test_fe_values = [&](const auto          &fe_values,
                                std::vector<double> &solution_values_out,
                                Vector<double>      &dst) {
        if (fe_values)
          {
            std::vector<Tensor<1, dim>> solution_gradients(
              fe_values->n_quadrature_points);
            std::vector<double> solution_values(fe_values->n_quadrature_points);

            for (const auto q : fe_values->quadrature_point_indices())
              {
                double         values = 0.;
                Tensor<1, dim> gradients;

                for (const auto i : fe_values->dof_indices())
                  {
                    gradients +=
                      solution_values_in[i] * fe_values->shape_grad(i, q);
                    values +=
                      solution_values_in[i] * fe_values->shape_value(i, q);
                  }
                solution_gradients[q] = gradients * fe_values->JxW(q);
                solution_values[q]    = values * fe_values->JxW(q);
              }

            for (const auto i : fe_values->dof_indices())
              {
                double sum_gradients = 0.;
                double sum_values    = 0.;
                for (const auto q : fe_values->quadrature_point_indices())
                  {
                    sum_gradients +=
                      solution_gradients[q] * fe_values->shape_grad(i, q);
                    sum_values +=
                      solution_values[q] * fe_values->shape_value(i, q);
                  }

                solution_values_cell_out[i] = sum_gradients + sum_values;
              }

            cell->distribute_local_to_global(
              Vector<double>(solution_values_cell_out.begin(),
                             solution_values_cell_out.end()),
              dst);
          }
      };

      test_fe_values(inside_fe_values, solution_values_cell_out, dst_cell_2);
      test_fe_values(surface_fe_values,
                     solution_values_surface_out,
                     dst_surface_2);

      for (const auto f : cell->face_indices())
        {
          if (cell->at_boundary(f))
            continue;

          non_matching_fe_interface_values.reinit(
            cell,
            f,
            numbers::invalid_unsigned_int,
            cell->neighbor(f),
            cell->neighbor_of_neighbor(f),
            numbers::invalid_unsigned_int);

          const auto &fe_interface_values =
            non_matching_fe_interface_values.get_inside_fe_values();

          if (fe_interface_values)
            {
              const auto &fe_face_values_m =
                fe_interface_values->get_fe_face_values(0);

              const auto &fe_face_values_p =
                fe_interface_values->get_fe_face_values(1);

              cell->neighbor(f)->get_dof_values(
                src,
                solution_values_neighbor_in.begin(),
                solution_values_neighbor_in.end());

              std::vector<Tensor<1, dim>> solution_gradients(
                fe_face_values_m.n_quadrature_points);
              std::vector<double> solution_values(
                fe_face_values_m.n_quadrature_points);

              for (const auto q : fe_face_values_m.quadrature_point_indices())
                {
                  double         values = 0.;
                  Tensor<1, dim> gradients;

                  for (const auto i : fe_face_values_m.dof_indices())
                    {
                      gradients += solution_values_in[i] *
                                     fe_face_values_m.shape_grad(i, q) -
                                   solution_values_neighbor_in[i] *
                                     fe_face_values_p.shape_grad(i, q);
                      values += solution_values_in[i] *
                                  fe_face_values_m.shape_value(i, q) -
                                solution_values_neighbor_in[i] *
                                  fe_face_values_p.shape_value(i, q);
                    }
                  solution_gradients[q] = gradients * fe_face_values_m.JxW(q);
                  solution_values[q]    = values * fe_face_values_m.JxW(q);
                }

              for (const auto i : fe_face_values_m.dof_indices())
                {
                  double sum_gradients = 0.;
                  double sum_values    = 0.;
                  for (const auto q :
                       fe_face_values_m.quadrature_point_indices())
                    {
                      sum_gradients += solution_gradients[q] *
                                       fe_face_values_m.shape_grad(i, q);
                      sum_values +=
                        solution_values[q] * fe_face_values_m.shape_value(i, q);
                    }

                  solution_values_faces_out[i] = sum_gradients + sum_values;
                }

              cell->distribute_local_to_global(
                Vector<double>(solution_values_faces_out.begin(),
                               solution_values_faces_out.end()),
                dst_faces_2);
            }
        }
    }

  deallog << "check difference l2 norm cell: "
          << dst_cell.l2_norm() - dst_cell_2.l2_norm() << std::endl;

  deallog << "check difference l2 norm surface: "
          << dst_surface.l2_norm() - dst_surface_2.l2_norm() << std::endl;

  deallog << "check difference l2 norm faces: "
          << dst_faces.l2_norm() - dst_faces_2.l2_norm() << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();

  test<2>(true);
  deallog << std::endl;
  test<2>(false);
  deallog << std::endl;
  test<3>(true);
  deallog << std::endl;
  test<3>(false);
}
